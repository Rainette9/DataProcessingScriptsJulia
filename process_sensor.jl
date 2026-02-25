# ── Sensor metadata ──────────────────────────────────────────────────────────

"""Column renames from raw TOA5 names to Peddy convention."""
const COLUMN_RENAME = Dict(
    :LI_H2Om   => :H2O,
    :LI_Pres   => :P,
    :diag_csat => :diag_sonic,
    :LI_diag   => :diag_gas,
)

"""
Per-sensor configuration: measurement height, slow-data column names,
whether a LI-7500 gas analyzer is present, and the column-name suffix
used in the raw fast-data files (empty string for SFC).
"""
const SENSORS = Dict(
    "SFC"    => (height=2,  ta_col=:TA,           rh_col=:RH,          has_licor=true,  suffix=""),
    "LOWER"  => (height=16, ta_col=:Temp_16m_Avg, rh_col=:RH_16m_Avg, has_licor=true,  suffix="_16m"),
    "UPPER"  => (height=26, ta_col=:Temp_26m_Avg, rh_col=:RH_26m_Avg, has_licor=true,  suffix="_26m"),
    "BOTTOM" => (height=5,  ta_col=:Temp_5m_Avg,  rh_col=:RH_5m_Avg,  has_licor=false, suffix="_5m"),
)

# ── Helper functions ─────────────────────────────────────────────────────────

"""
    build_hf_dimarray(fast_data::DataFrame, meta) -> DimArray

Build a high-frequency DimArray from raw fast data, selecting columns based on
whether the sensor has a LICOR and renaming via `COLUMN_RENAME`.
Handles per-sensor column suffixes (e.g. `Ux_26m` → `Ux` for UPPER).
"""
function build_hf_dimarray(fast_data::DataFrame, meta)
    if meta.has_licor
        base_cols = [:Ux, :Uy, :Uz, :Ts, :LI_H2Om, :LI_Pres, :diag_csat, :LI_diag]
    else
        base_cols = [:Ux, :Uy, :Uz, :Ts, :diag_csat]
    end

    # Append sensor suffix to get the actual DataFrame column names
    s = meta.suffix
    raw_cols = [Symbol(string(c) * s) for c in base_cols]

    # Peddy variable names: apply COLUMN_RENAME to the base (unsuffixed) names
    hf_vars = [get(COLUMN_RENAME, c, c) for c in base_cols]
    hf_matrix = hcat([Float64.(fast_data[!, c]) for c in raw_cols]...)

    return DimArray(hf_matrix, (Ti(fast_data.TIMESTAMP), Var(hf_vars)))
end

"""
    build_lf_dimarray(slow_data::DataFrame, meta, date::Date) -> Union{DimArray, Nothing}

Filter slow data to a single day and build a low-frequency DimArray with TA and
RH columns. RH is divided by 100 to convert from percentage to fraction.
Returns `nothing` if no data is available for the given day.
"""
function build_lf_dimarray(slow_data::DataFrame, meta, date::Date)
    ta_name = string(meta.ta_col)
    rh_name = string(meta.rh_col)

    if !(ta_name in names(slow_data) && rh_name in names(slow_data))
        return nothing
    end

    day_start = DateTime(date)
    day_end   = DateTime(date + Day(1))
    day_slow  = filter(r -> r.TIMESTAMP >= day_start && r.TIMESTAMP < day_end, slow_data)

    if nrow(day_slow) == 0
        return nothing
    end

    return DimArray(
        hcat(
            Float64.(coalesce.(day_slow[!, ta_name], NaN)),
            Float64.(coalesce.(day_slow[!, rh_name], NaN)) ./ 100  # fraction 0–1
        ),
        (Ti(day_slow.TIMESTAMP), Var([:TA, :RH]))
    )
end

"""
    processed_dates(output_dir::String, sensor::String) -> Set{Date}

Scan the output directory for already-processed hourly CSV files and return the
set of dates that have at least one output file.
"""
function processed_dates(output_dir::String, sensor::String)
    sensor_dir = joinpath(output_dir, sensor)
    if !isdir(sensor_dir)
        return Set{Date}()
    end

    pattern = Regex("$(sensor)_Fastdata_proc_(\\d{4}-\\d{2}-\\d{2})_\\d{4}\\.dat\$")
    found = Set{Date}()

    for (root, _, files) in walkdir(sensor_dir)
        for f in files
            m = match(pattern, f)
            if m !== nothing
                push!(found, Date(m.captures[1], "yyyy-mm-dd"))
            end
        end
    end

    return found
end

"""
    make_pipeline(meta, cal_coeffs, has_slow; spike_threshold, window_minutes,
                  block_duration_min, max_gap_minutes, dt_ms, logger)

Construct the Peddy sensor, variable groups, and EddyPipeline.
Returns `(pipeline, output)`.
"""
function make_pipeline(meta, cal_coeffs, has_slow;
                       spike_threshold, window_minutes, block_duration_min,
                       max_gap_minutes, dt_ms, logger)
    peddy_sensor = meta.has_licor ? LICOR(calibration_coefficients=cal_coeffs) : CSAT3()

    has_gas_analyzer = meta.has_licor && cal_coeffs !== nothing && has_slow

    var_groups = if meta.has_licor
        [
            VariableGroup("Wind",        [:Ux, :Uy, :Uz], spike_threshold=spike_threshold),
            VariableGroup("Temperature", [:Ts],            spike_threshold=spike_threshold),
            VariableGroup("Gas",         [:H2O],           spike_threshold=spike_threshold),
        ]
    else
        [
            VariableGroup("Wind",        [:Ux, :Uy, :Uz], spike_threshold=spike_threshold),
            VariableGroup("Temperature", [:Ts],            spike_threshold=spike_threshold),
        ]
    end

    output = MemoryOutput()
    pipeline = EddyPipeline(
        sensor          = peddy_sensor,
        quality_control = PhysicsBoundsCheck(
            Ux = Peddy.Limit(-80.0, 80.0),
            Uy = Peddy.Limit(-80.0, 80.0),
            Uz = Peddy.Limit(-80.0, 80.0),
            Ts = Peddy.Limit(-50.0, 20.0),
        ),
        gas_analyzer    = has_gas_analyzer ? H2OCalibration() : nothing,
        despiking       = SimpleSigmundDespiking(
            window_minutes  = window_minutes,
            variable_groups = var_groups,
        ),
        make_continuous = MakeContinuous(step_size_ms=dt_ms, max_gap_minutes=max_gap_minutes),
        gap_filling     = GeneralInterpolation(),
        double_rotation = WindDoubleRotation(block_duration_minutes=block_duration_min),
        output          = output,
        logger          = logger,
    )

    return pipeline, output
end

# ── Main processing function ─────────────────────────────────────────────────

"""
    process_sensor(; input_base, sensor, processed_output, year=2025,
                    dates=nothing, resume=true, config=Dict(),
                    slow_data=DataFrame(), location=LocationMetadata(...))

Batch-process a full year of eddy covariance data for one sensor, day by day.

Reads raw TOA5 fast data, runs the Peddy.jl pipeline, and writes hourly
processed fast-data CSVs to `processed_output/{sensor}/{YYYYmm}/`.

Slow data must be read, cleaned, and saved separately (see the notebook's
slow-data cell) and passed in via the `slow_data` keyword argument.
When `slow_data` is empty, the gas analyzer correction is skipped.

When `resume=true` (default), dates that already have output files are skipped.
Set `resume=false` to reprocess all dates from scratch (existing files are
overwritten).

# Config Dict Keys (all optional)
| Key                      | Type    | Default                              |
|--------------------------|---------|--------------------------------------|
| `:sensor_type`           | Symbol  | `:licor` for SFC/LOWER/UPPER, `:csat3` for BOTTOM |
| `:spike_threshold`       | Float64 | 6.0                                  |
| `:window_minutes`        | Float64 | 5.0                                  |
| `:block_duration_minutes`| Float64 | 30.0                                 |
| `:max_gap_minutes`       | Float64 | 5.0                                  |
"""
function process_sensor(;
    input_base::String,
    sensor::String,
    processed_output::String,
    year::Int = 2025,
    dates = nothing,
    resume::Bool = true,
    config::Dict = Dict(),
    slow_data::DataFrame = DataFrame(),
    location = LocationMetadata(
        latitude          = -71.949608,
        longitude         = 23.0,
        elevation         = 1382.0,
        instrument_height = Float64(SENSORS[sensor].height),
    )
)
    # --- Look up sensor metadata ---
    meta = SENSORS[sensor]

    spike_threshold    = get(config, :spike_threshold, 6.0)
    window_minutes     = get(config, :window_minutes, 5.0)
    block_duration_min = get(config, :block_duration_minutes, 30.0)
    max_gap_minutes    = get(config, :max_gap_minutes, 5.0)

    sensor_type = meta.has_licor ? :licor : :csat3

    println("=" ^ 60)
    println("Processing sensor: $sensor  (type: $sensor_type)")
    println("  TA column: $(meta.ta_col)  |  RH column: $(meta.rh_col)")
    has_slow = nrow(slow_data) > 0
    if has_slow
        println("  Slow data: $(nrow(slow_data)) records")
    else
        println("  Slow data: none provided — gas analyzer correction will be skipped")
    end
    println("=" ^ 60)

    # --- Determine dates to process ---
    if isnothing(dates)
        println("\n--- Listing available dates ---")
        dates = list_available_dates(base_path=input_base, sensor=sensor)
    end

    if isempty(dates)
        @warn "No dates to process for sensor $sensor"
        return nothing
    end

    # --- Skip already-processed dates (when resuming) ---
    if resume
        done = processed_dates(processed_output, sensor)
        remaining = filter(d -> d ∉ done, dates)
        n_skipped = length(dates) - length(remaining)

        if n_skipped > 0
            println("\n  Already processed: $n_skipped dates (skipped)")
        end

        if isempty(remaining)
            println("\nAll dates already processed — nothing to do.")
            return nothing
        end
    else
        remaining = dates
        println("\n  Reprocessing all dates from scratch")
    end

    println("  Dates to process: $(length(remaining))  " *
            "($(first(remaining)) to $(last(remaining)))")

    # --- Calibration coefficients (once) ---
    cal_coeffs = meta.has_licor ?
        Peddy.default_calibration_coefficients(sensor, year; number_type=Float64) : nothing

    # --- Output directory ---
    fast_out_base = joinpath(processed_output, sensor)
    mkpath(fast_out_base)

    # --- Processing logger ---
    logger = ProcessingLogger()
    log_metadata!(logger, "sensor", sensor)
    log_metadata!(logger, "sensor_type", string(sensor_type))
    log_metadata!(logger, "year", string(year))

    has_gas_analyzer = meta.has_licor && cal_coeffs !== nothing && has_slow
    log_metadata!(logger, "gas_analyzer", has_gas_analyzer ? "H2OCalibration" : "disabled")
    if cal_coeffs !== nothing
        log_metadata!(logger, "calibration_A", string(cal_coeffs.A))
        log_metadata!(logger, "calibration_B", string(cal_coeffs.B))
        log_metadata!(logger, "calibration_C", string(cal_coeffs.C))
        log_metadata!(logger, "calibration_H2O_Zero", string(cal_coeffs.H2O_Zero))
        log_metadata!(logger, "calibration_H20_Span", string(cal_coeffs.H20_Span))
    end

    current_month = nothing
    log_file = ""

    # --- Day-by-day processing loop ---
    total_files = 0
    failed_days = Date[]

    for (i, date) in enumerate(remaining)
        println("\n--- [$i/$(length(remaining))] Processing $date ---")

        try
            # Read fast data for this day
            fast_data = read_fast_data(base_path=input_base, sensor=sensor, dates=[date])

            if nrow(fast_data) == 0
                @warn "No fast data for $date — skipping"
                continue
            end

            unique!(fast_data, :TIMESTAMP)
            sort!(fast_data, :TIMESTAMP)

            # Detect sampling rate
            dt_ms = round(Int, median([Dates.value(d) for d in diff(fast_data.TIMESTAMP)]))
            println("  Records: $(nrow(fast_data))  Sampling: $(dt_ms) ms ($(round(1000.0/dt_ms, digits=1)) Hz)")

            # Build DimArrays
            high_frequency_data = build_hf_dimarray(fast_data, meta)
            low_frequency_data  = has_slow ? build_lf_dimarray(slow_data, meta, date) : nothing

            if low_frequency_data === nothing && has_slow
                println("  No slow data for $date — skipping gas analyzer correction")
            end

            # Build pipeline (disable gas analyzer when no slow data for this day)
            day_has_slow = low_frequency_data !== nothing
            pipeline, output = make_pipeline(meta, cal_coeffs, day_has_slow;
                spike_threshold, window_minutes, block_duration_min,
                max_gap_minutes, dt_ms, logger)

            process!(pipeline, high_frequency_data, low_frequency_data)

            # Extract results and write hourly dat files
            hf_result, _ = Peddy.get_results(output)
            n_written = _write_hourly_files(hf_result, fast_out_base, sensor, dt_ms)
            total_files += n_written
            println("  Wrote $n_written hourly dat files")

            # Write/append to monthly processing log
            month_str = Dates.format(date, "yyyymm")
            if month_str != current_month
                current_month = month_str
                month_dir = joinpath(fast_out_base, month_str)
                mkpath(month_dir)
                log_file = joinpath(month_dir, "$(sensor)_processing_log_$(month_str).csv")
                write_processing_log(logger, log_file)
            else
                write_processing_log(logger, log_file; append=true)
            end
            reset!(logger)


        catch e
            push!(failed_days, date)
            @warn "Failed to process $date" exception=(e, catch_backtrace())
        end

        # Free memory
        GC.gc()
    end

    # --- Summary ---
    println("\n" * "=" ^ 60)
    println("Finished processing sensor: $sensor")
    println("  Total hourly files written: $total_files")
    println("  Failed days: $(length(failed_days))")
    if !isempty(failed_days)
        println("  Failed dates: ", join(string.(failed_days), ", "))
    end
    println("=" ^ 60)

    return nothing
end

# ── Output writer ────────────────────────────────────────────────────────────

"""
    _write_hourly_files(hf_result::DimArray, out_base::String, sensor::String, dt_ms::Int) -> Int

Write one `.dat` file per hour covered by the processed DimArray. Every file
contains a complete hour of timestamps (from XX:00:00.000 to XX:59:59.xxx) at
the sampling interval `dt_ms`, so all files have exactly `div(3_600_000, dt_ms)`
rows. Timestamps without data are filled with NaN.

Returns the number of files written.
"""
function _write_hourly_files(hf_result::DimArray, out_base::String, sensor::String, dt_ms::Int)
    hf_times = collect(dims(hf_result, Ti))
    hf_vars  = collect(dims(hf_result, Var))

    if isempty(hf_times)
        return 0
    end

    rows_per_hour = div(3_600_000, dt_ms)
    step = Millisecond(dt_ms)

    # Extract the underlying data matrix once (rows = time, cols = vars)
    raw_matrix = parent(hf_result)

    # Build a lookup: timestamp → row index in the raw matrix
    time_index = Dict{DateTime, Int}()
    for (i, t) in enumerate(hf_times)
        time_index[t] = i
    end

    # Determine which full hours are covered
    first_hour = floor(hf_times[1],  Hour)
    last_hour  = floor(hf_times[end], Hour)

    file_count = 0
    t_hour = first_hour
    n_vars = length(hf_vars)

    while t_hour <= last_hour
        # Build the complete time axis for this hour
        full_times = [t_hour + (k - 1) * step for k in 1:rows_per_hour]

        # Pre-fill with NaN
        matrix = fill(NaN, rows_per_hour, n_vars)

        # Copy available data into the matrix
        has_data = false
        for (row, t) in enumerate(full_times)
            idx = get(time_index, t, 0)
            if idx > 0
                @views matrix[row, :] .= raw_matrix[idx, :]
                has_data = true
            end
        end

        # Skip hours with no data at all
        if !has_data
            t_hour += Hour(1)
            continue
        end

        # Build DataFrame
        df = DataFrame(:timestamp => full_times)
        for (col, v) in enumerate(hf_vars)
            df[!, v] = @view matrix[:, col]
        end

        # Create monthly subfolder and write
        month_dir = joinpath(out_base, Dates.format(Date(t_hour), "yyyymm"))
        mkpath(month_dir)

        fname = joinpath(month_dir,
            "$(sensor)_Fastdata_proc_$(Dates.format(t_hour, "yyyy-mm-dd_HHMM")).dat")
        CSV.write(fname, df; missingstring="NaN")
        file_count += 1

        t_hour += Hour(1)
    end

    return file_count
end
