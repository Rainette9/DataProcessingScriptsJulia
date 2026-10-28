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
whether a LI-7500 gas analyzer is present, the column-name suffix
used in the raw fast-data files (empty string for SFC), the offset
added to the raw LI-7500 pressure (the SFC logger stores kPa - 80), and the
number of samples on either side of a CSAT warning that are also removed.
The margin is 0 at SFC: there a warning marks a single sample blocked by a snow
particle and its neighbours are as good as any other sample. On the towers the
samples next to a warning are clearly disturbed (sonic faults), so they get 2.
"""
const SENSORS = Dict(
    "SFC"    => (height=2,  ta_col=:TA,           rh_col=:RH,          has_licor=true,  suffix="",     pres_offset=80.0, diag_margin=0),
    "LOWER"  => (height=16, ta_col=:Temp_16m_Avg, rh_col=:RH_16m_Avg, has_licor=true,  suffix="_16m", pres_offset=0.0,  diag_margin=2),
    "UPPER"  => (height=26, ta_col=:Temp_26m_Avg, rh_col=:RH_26m_Avg, has_licor=true,  suffix="_26m", pres_offset=0.0,  diag_margin=2),
    "BOTTOM" => (height=5,  ta_col=:Temp_5m_Avg,  rh_col=:RH_5m_Avg,  has_licor=false, suffix="_5m",  pres_offset=0.0,  diag_margin=2),
)

"""Physically plausible H2O molar density (mmol m⁻³) after the gas analyzer correction."""
const H2O_BOUNDS = (0.0, 400.0)

"""Reference (slow data) gaps up to this many minutes are interpolated for the H2O correction."""
const MAX_REF_GAP_MIN = 10

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
    # Raw files where a column is only "NAN" are read as String columns
    to_float(x) = ismissing(x) ? NaN :
                  x isa AbstractString ? something(tryparse(Float64, x), NaN) : Float64(x)
    hf_matrix = hcat([map(to_float, fast_data[!, c]) for c in raw_cols]...)

    # Absolute pressure in kPa (needed by the H2O correction)
    ip = findfirst(==(:P), hf_vars)
    if ip !== nothing && meta.pres_offset != 0
        hf_matrix[:, ip] .+= meta.pres_offset
    end

    return DimArray(hf_matrix, (Ti(fast_data.TIMESTAMP), Var(hf_vars)))
end

"""Linear interpolation over NaN runs of at most `maxgap` points; longer runs stay NaN."""
function fill_short_gaps!(x::Vector{Float64}, maxgap::Int)
    n = length(x)
    i = 1
    while i <= n
        if isnan(x[i])
            j = i
            while j < n && isnan(x[j + 1])
                j += 1
            end
            if i > 1 && j < n && (j - i + 1) <= maxgap
                for k in i:j
                    w = (k - i + 1) / (j - i + 2)
                    x[k] = (1 - w) * x[i - 1] + w * x[j + 1]
                end
            end
            i = j + 1
        else
            i += 1
        end
    end
    return x
end

"""
    _lf_grid(times, ta_vals, rh_vals, date::Date) -> Union{DimArray, Nothing}

Place TA (°C) / RH (%) samples on a regular 1-min grid from `date - 1h` to
`date + 1d + 1h`. Minutes without a sample are NaN, so the H2O correction gives
NaN there instead of using a reference minute hours away (Peddy maps each fast
sample to the nearest slow timestamp without a distance limit). Gaps up to
`MAX_REF_GAP_MIN` minutes are interpolated. Logger zeros (TA and RH both exactly
0) are treated as missing. RH is divided by 100 to the fraction Peddy expects.
Returns `nothing` if the day itself has no valid TA/RH minute.
"""
function _lf_grid(times, ta_vals, rh_vals, date::Date)
    grid = collect(DateTime(date) - Hour(1):Minute(1):DateTime(date + Day(1)) + Hour(1))

    lookup = Dict{DateTime, Tuple{Float64, Float64}}()
    for (t, a, b) in zip(times, ta_vals, rh_vals)
        (grid[1] <= t <= grid[end]) || continue
        ta = ismissing(a) ? NaN : Float64(a)
        rh = ismissing(b) ? NaN : Float64(b)
        if ta == 0 && rh == 0
            continue  # logger zeros
        end
        lookup[floor(t, Minute)] = (ta, rh / 100)  # fraction 0–1
    end

    ta = [get(lookup, t, (NaN, NaN))[1] for t in grid]
    rh = [get(lookup, t, (NaN, NaN))[2] for t in grid]

    day = DateTime(date) .<= grid .< DateTime(date + Day(1))
    if !any(.!isnan.(rh[day]) .& .!isnan.(ta[day]))
        return nothing
    end

    fill_short_gaps!(ta, MAX_REF_GAP_MIN)
    fill_short_gaps!(rh, MAX_REF_GAP_MIN)

    return DimArray(hcat(ta, rh), (Ti(grid), Var([:TA, :RH])))
end

"""
    build_lf_dimarray(slow_data::DataFrame, meta, date::Date) -> Union{DimArray, Nothing}

Own-probe path (SFC, or towers without a reference chain): reference from the
sensor's `meta.ta_col` / `meta.rh_col` (RH in percent). See `_lf_grid`.
"""
function build_lf_dimarray(slow_data::DataFrame, meta, date::Date)
    ta_name = string(meta.ta_col)
    rh_name = string(meta.rh_col)
    (ta_name in names(slow_data) && rh_name in names(slow_data)) || return nothing
    return _lf_grid(slow_data.TIMESTAMP, slow_data[!, ta_name], slow_data[!, rh_name], date)
end

"""
    build_lf_dimarray(ref::DataFrame, date::Date) -> Union{DimArray, Nothing}

Reference-chain path (LOWER/UPPER): `ref` has `TIMESTAMP`, `TA_ref`, `RH_ref`
(RH in percent) and `REF_SRC`, see `build_reference_chain` in `read_data,jl`.
"""
build_lf_dimarray(ref::DataFrame, date::Date) = _lf_grid(ref.TIMESTAMP, ref.TA_ref, ref.RH_ref, date)

"""
    _ref_source_counts(ref::DataFrame, date::Date) -> NamedTuple

Minutes of the day by reference source (chain entries 1–4, 9 = interpolated, 0 = none).
"""
function _ref_source_counts(ref::DataFrame, date::Date)
    d0 = DateTime(date); d1 = DateTime(date + Day(1))
    src = ref.REF_SRC[(ref.TIMESTAMP .>= d0) .& (ref.TIMESTAMP .< d1)]
    c(k) = count(==(k), src)
    return (src1=c(1), src2=c(2), src3=c(3), src4=c(4), interp=c(9), none=c(0) + (1440 - length(src)))
end

"""View of the H2O column of a pipeline result, whatever its (Ti, Var) orientation; `nothing` if absent."""
function _h2o_column(hf_result::DimArray)
    ih = findfirst(==(:H2O), collect(dims(hf_result, Var)))
    ih === nothing && return nothing
    M = parent(hf_result)
    return dimnum(hf_result, Var) == 2 ? (@view M[:, ih]) : (@view M[ih, :])
end

"""Set H2O outside `H2O_BOUNDS` to NaN (in place). Returns the number of values removed."""
function apply_h2o_bounds!(hf_result::DimArray)
    col = _h2o_column(hf_result)
    col === nothing && return 0
    bad = .!isnan.(col) .& ((col .< H2O_BOUNDS[1]) .| (col .> H2O_BOUNDS[2]))
    col[bad] .= NaN
    return sum(bad)
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
                  block_duration_min, max_gap_minutes, max_gap_samples, diag_margin,
                  ts_max_deviation, dt_ms, logger, gap_filling=true)

Construct the Peddy sensor, variable groups, and EddyPipeline.
Returns `(pipeline, output)`.
"""
function make_pipeline(meta, cal_coeffs, has_slow;
                       spike_threshold, window_minutes, block_duration_min,
                       max_gap_minutes, max_gap_samples, diag_margin, ts_max_deviation,
                       dt_ms, logger, gap_filling=true)
    peddy_sensor = meta.has_licor ?
        LICOR(calibration_coefficients=cal_coeffs, diag_margin=diag_margin) :
        CSAT3(diag_margin=diag_margin)

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
            max_ts_deviation = ts_max_deviation,
            couple_sonic     = true,  # a despiked Ts or wind value removes the whole sonic sample
        ),
        make_continuous = MakeContinuous(step_size_ms=dt_ms, max_gap_minutes=max_gap_minutes),
        gap_filling     = gap_filling ? GeneralInterpolation(max_gap_size=max_gap_samples) : nothing,
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
                    slow_data=DataFrame(), ref_data=nothing, log_tag="",
                    location=LocationMetadata(...))

Batch-process a full year of eddy covariance data for one sensor, day by day.

Reads raw TOA5 fast data, runs the Peddy.jl pipeline, and writes hourly
processed fast-data CSVs to `processed_output/{sensor}/{YYYYmm}/`.

Slow data must be read, cleaned, and saved separately and passed in via the
`slow_data` keyword argument. When `slow_data` is empty, the gas analyzer
correction is skipped.

`ref_data` (optional) is a minute-by-minute TA/RH reference built with
`build_reference_chain` (columns `TIMESTAMP`, `TA_ref`, `RH_ref`, `REF_SRC`).
When given, it replaces `slow_data` as the source of the H2O correction
reference, and a per-day source summary is written to
`{processed_output}/{sensor}/{YYYYmm}/{sensor}_h2o_ref_{YYYYmm}.csv`. Days
without any reference minute get NaN H2O (the raw LI-7500 H2O must not pass
through uncorrected).

After the pipeline, H2O outside `H2O_BOUNDS` (0–400 mmol m⁻³) is set to NaN.

When `resume=true` (default), dates that already have output files are skipped
and the monthly processing log is appended to. Set `resume=false` to reprocess
all dates from scratch (existing files and logs are overwritten).

`log_tag` is appended to the monthly processing log names (e.g. `_patch`) so a
partial run does not overwrite the logs of a full run.

# Config Dict Keys (all optional)
| Key                      | Type    | Default                              |
|--------------------------|---------|--------------------------------------|
| `:sensor_type`           | Symbol  | `:licor` for SFC/LOWER/UPPER, `:csat3` for BOTTOM |
| `:spike_threshold`       | Float64 | 6.0                                  |
| `:window_minutes`        | Float64 | 5.0                                  |
| `:block_duration_minutes`| Float64 | 30.0                                 |
| `:max_gap_minutes`       | Float64 | 5.0                                  |
| `:max_gap_samples`       | Int     | 5 (longest gap filled by interpolation, in samples) |
| `:gap_filling`           | Bool    | true (false: removed samples stay NaN) |
| `:diag_margin`           | Int     | per sensor from `SENSORS` (SFC 0, towers 2; samples either side of a CSAT warning also set to NaN) |
| `:ts_max_deviation`      | Float64 | 10.0 (K; |Ts - 5-min rolling median| above this removes the sonic sample) |
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
    ref_data::Union{Nothing, DataFrame} = nothing,
    log_tag::String = "",
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
    max_gap_samples    = get(config, :max_gap_samples, 5)
    gap_filling        = get(config, :gap_filling, true)
    diag_margin        = get(config, :diag_margin, meta.diag_margin)
    ts_max_deviation   = get(config, :ts_max_deviation, 10.0)

    sensor_type = meta.has_licor ? :licor : :csat3

    println("=" ^ 60)
    println("Processing sensor: $sensor  (type: $sensor_type)")
    println("  TA column: $(meta.ta_col)  |  RH column: $(meta.rh_col)")
    use_ref  = ref_data !== nothing && nrow(ref_data) > 0
    has_slow = use_ref || nrow(slow_data) > 0
    if use_ref
        println("  H2O reference: chain, $(nrow(ref_data)) minutes " *
                "($(ref_data.TIMESTAMP[1]) → $(ref_data.TIMESTAMP[end]))")
    elseif has_slow
        println("  Slow data: $(nrow(slow_data)) records")
    else
        println("  Slow data: none provided — gas analyzer correction will be skipped")
    end
    println("  diag_margin=$diag_margin  ts_max_deviation=$ts_max_deviation K  " *
            "gap_filling=$gap_filling (max $max_gap_samples samples)")
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
    log_metadata!(logger, "interpolation_max_gap_samples", string(max_gap_samples))
    log_metadata!(logger, "gap_filling", string(gap_filling))
    log_metadata!(logger, "sonic_diag_margin_samples", string(diag_margin))
    log_metadata!(logger, "ts_max_deviation_K", string(ts_max_deviation))
    log_metadata!(logger, "despiking_couple_sonic", "true")
    log_metadata!(logger, "pressure_offset_kPa", string(meta.pres_offset))
    log_metadata!(logger, "h2o_reference", use_ref ? "chain" : (has_slow ? "own_slow_data" : "none"))

    has_gas_analyzer = meta.has_licor && cal_coeffs !== nothing && has_slow
    log_metadata!(logger, "gas_analyzer", has_gas_analyzer ? "H2OCalibration" : "disabled")
    if cal_coeffs !== nothing
        log_metadata!(logger, "calibration_A", string(cal_coeffs.A))
        log_metadata!(logger, "calibration_B", string(cal_coeffs.B))
        log_metadata!(logger, "calibration_C", string(cal_coeffs.C))
        log_metadata!(logger, "calibration_H2O_Zero", string(cal_coeffs.H2O_Zero))
        log_metadata!(logger, "calibration_H20_Span", string(cal_coeffs.H20_Span))
    end

    # Months whose logs this process has already opened. On the first touch of a
    # month the processing log is overwritten unless resuming onto an existing
    # file (the daily cron runs in resume mode and must not wipe earlier days).
    months_started = Set{String}()

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

            # Keep only this calendar day: the raw daily files overlap (each spans ~25 h and
            # the reader also loads the neighbouring file), so without this filter every run
            # reprocesses ~2 days and overwrites the neighbouring day's hourly output.
            n_all = nrow(fast_data)
            filter!(r -> Date(r.TIMESTAMP) == date, fast_data)
            n_all > nrow(fast_data) && println("  Kept $(nrow(fast_data)) of $n_all records belonging to $date")
            if nrow(fast_data) == 0
                @warn "No fast data for $date after trimming to the day — skipping"
                continue
            end

            # Detect sampling rate
            dt_ms = round(Int, median([Dates.value(d) for d in diff(fast_data.TIMESTAMP)]))
            println("  Records: $(nrow(fast_data))  Sampling: $(dt_ms) ms ($(round(1000.0/dt_ms, digits=1)) Hz)")

            # Build DimArrays
            high_frequency_data = build_hf_dimarray(fast_data, meta)
            low_frequency_data  = use_ref  ? build_lf_dimarray(ref_data, date) :
                                  has_slow ? build_lf_dimarray(slow_data, meta, date) : nothing
            rc = use_ref ? _ref_source_counts(ref_data, date) : nothing
            if use_ref
                println("  Reference minutes: src1=$(rc.src1) src2=$(rc.src2) src3=$(rc.src3) " *
                        "src4=$(rc.src4) interp=$(rc.interp) none=$(rc.none)")
            end

            no_reference_day = low_frequency_data === nothing && has_slow
            if no_reference_day
                println("  No slow data for $date — gas analyzer correction skipped" *
                        (use_ref ? "; H2O will be set to NaN (no reference available)" : ""))
            end

            # Build pipeline (disable gas analyzer when no slow data for this day)
            day_has_slow = low_frequency_data !== nothing
            pipeline, output = make_pipeline(meta, cal_coeffs, day_has_slow;
                spike_threshold, window_minutes, block_duration_min,
                max_gap_minutes, max_gap_samples, diag_margin, ts_max_deviation,
                dt_ms, logger, gap_filling)

            process!(pipeline, high_frequency_data, low_frequency_data)

            # Extract results and write hourly dat files
            hf_result, _ = Peddy.get_results(output)

            # Remove physically impossible H2O (mmol m-3) after the gas analyzer correction
            if meta.has_licor
                n_bad_h2o = apply_h2o_bounds!(hf_result)
                n_bad_h2o > 0 && println("  Removed $n_bad_h2o H2O values outside $(H2O_BOUNDS) mmol m-3")
                log_event!(logger, :postprocess, :h2o_bounds; variable=:H2O, count=n_bad_h2o,
                           note="samples outside $(H2O_BOUNDS[1])-$(H2O_BOUNDS[2]) mmol/m3 set to NaN")
                if use_ref && no_reference_day
                    # Reference chain in use but no reference minute for this day: the raw
                    # LI-7500 H2O must not pass through uncorrected.
                    col = _h2o_column(hf_result)
                    col === nothing || (col .= NaN)
                    println("  H2O set to NaN for the whole day (no reference)")
                    log_event!(logger, :postprocess, :h2o_no_reference; variable=:H2O,
                               note="no TA/RH reference for this day; H2O set to NaN")
                end
            end

            n_written = _write_hourly_files(hf_result, fast_out_base, sensor, dt_ms)
            total_files += n_written
            println("  Wrote $n_written hourly dat files")

            # Write/append to monthly processing log (+ reference-source summary)
            month_str = Dates.format(date, "yyyymm")
            month_dir = joinpath(fast_out_base, month_str)
            mkpath(month_dir)
            log_file = joinpath(month_dir, "$(sensor)_processing_log_$(month_str)$(log_tag).csv")
            ref_file = joinpath(month_dir, "$(sensor)_h2o_ref_$(month_str)$(log_tag).csv")
            first_touch = !(month_str in months_started)
            push!(months_started, month_str)
            append_log = !first_touch || (resume && isfile(log_file))
            write_processing_log(logger, log_file; append=append_log)
            if use_ref
                append_ref = !first_touch || (resume && isfile(ref_file))
                open(ref_file, append_ref ? "a" : "w") do io
                    append_ref || println(io, "date,src1,src2,src3,src4,interp,none")
                    println(io, "$date,$(rc.src1),$(rc.src2),$(rc.src3),$(rc.src4),$(rc.interp),$(rc.none)")
                end
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
