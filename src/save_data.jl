using CSV
using DataFrames
using Dates

# Unit lookup for cleaned slow data columns (matching original TOA5 format)
const SLOW_DATA_UNITS = Dict(
    "TIMESTAMP" => "TS",
    "RECORD"    => "RN",
    "WD1"       => "deg",
    "WD2"       => "deg",
    "TA"        => "degC",
    "RH"        => "%",
    "HS_Cor"    => "m",
    "SFTempK"   => "K",
    "SWdown1"   => "W/m2",
    "SWdown2"   => "W/m2",
    "SWup1"     => "W/m2",
    "SWup2"     => "W/m2",
    "LWdown1"   => "W/m2",
    "LWdown2"   => "W/m2",
    "LWup1"     => "W/m2",
    "LWup2"     => "W/m2",
    "SWdn"      => "W/m2",
    "SensorT"   => "degC",
    "PF_FC4"    => "g/m2/s",
    "WS_FC4"    => "m/s",
    "WS1_Avg"   => "m/s",
    "WS2_Avg"   => "m/s",
    "WS1_Max"   => "m/s",
    "WS2_Max"   => "m/s",
    "WS1_Std"   => "m/s",
    "WS2_Std"   => "m/s",
    "Temp_5m_Avg"  => "degC",
    "RH_5m_Avg"    => "%",
    "Temp_16m_Avg" => "degC",
    "RH_16m_Avg"   => "%",
    "Temp_26m_Avg" => "degC",
    "RH_26m_Avg"   => "%",
    "SBTempK"      => "K",
    "Ux_Avg"       => "m/s",
    "Uy_Avg"       => "m/s",
    "Uz_Avg"       => "m/s",
    "Ts_Avg"       => "degC",
    "LI_Pres_Avg"  => "kPa",
    # Tower sensor columns
    "VP_16m_Avg"   => "kPa",
    "AP_16m_Avg"   => "kPa",
    "Temp_10m_Avg" => "degC",
    "RH_10m_Avg"   => "%",
    "WS_16m_Avg"   => "m/s",
    "WS_16m_Max"   => "m/s",
    "WS_16m_Std"   => "m/s",
    "WD_16m"       => "deg",
    "WS_10m_Avg"   => "m/s",
    "WS_10m_Max"   => "m/s",
    "WS_10m_Std"   => "m/s",
    "WD_10m"       => "deg",
    "FluxMin_11m"  => "g/m2/s",
    "FluxMean_11m" => "g/m2/s",
    "FluxMax_11m"  => "g/m2/s",
    "FluxStd_11m"  => "g/m2/s",
    "FluxSum_11m"  => "g/m2/s",
    "WindMin_11m"  => "m/s",
    "WindMean_11m" => "m/s",
    "WindMax_11m"  => "m/s",
    "FC_Flux_11m"  => "g/m2",
    "FC_Wind_11m"  => "m/s",
    "FluxMin_4m"   => "g/m2/s",
    "FluxMean_4m"  => "g/m2/s",
    "FluxMax_4m"   => "g/m2/s",
    "FluxStd_4m"   => "g/m2/s",
    "FluxSum_4m"   => "g/m2/s",
    "WindMin_4m"   => "m/s",
    "WindMean_4m"  => "m/s",
    "WindMax_4m"   => "m/s",
    "FC_Flux_4m"   => "g/m2",
    "FC_Wind_4m"   => "m/s",
    "Incoming_SW_16m_Avg"  => "W/m2",
    "Outgoing_SW_16m_Avg"  => "W/m2",
    "Incoming_LW_16m_Avg"  => "W/m2",
    "Outgoing_LW_16m_Avg"  => "W/m2",
    "Ux_16m_Avg"   => "m/s",
    "Uy_16m_Avg"   => "m/s",
    "Uz_16m_Avg"   => "m/s",
    "Ts_16m_Avg"   => "degC",
    "LI_H2Om_16m_Avg" => "mmol/m3",
    "LI_Pres_16m_Avg" => "kPa",
    "VP_26m_Avg"   => "kPa",
    "AP_26m_Avg"   => "kPa",
    "WS_26m_Avg"   => "m/s",
    "WS_26m_Max"   => "m/s",
    "WS_26m_Std"   => "m/s",
    "WD_26m"       => "deg",
    "Incoming_SW_26m_Avg"  => "W/m2",
    "Outgoing_SW_26m_Avg"  => "W/m2",
    "Incoming_UW_26m_Avg"  => "W/m2",
    "Outgoing_UW_26m_Avg"  => "W/m2",
    "Ux_26m_Avg"   => "m/s",
    "Uy_26m_Avg"   => "m/s",
    "Uz_26m_Avg"   => "m/s",
    "Ts_26m_Avg"   => "degC",
    "LI_H2Om_26m_Avg" => "mmol/m3",
    "LI_Pres_26m_Avg" => "kPa",
    "WS_5m_Avg"    => "m/s",
    "WS_5m_Max"    => "m/s",
    "WS_5m_Std"    => "m/s",
    "Ux_5m_Avg"    => "m/s",
    "Uy_5m_Avg"    => "m/s",
    "Uz_5m_Avg"    => "m/s",
    "Ts_5m_Avg"    => "degC",
    "SWdn_5m"      => "W/m2",
    "CS320_T_5m"   => "degC",
    "wnd_spd_5m"   => "m/s",
)

# Processing type lookup (matching original TOA5 row 4)
const SLOW_DATA_PROCESSING = Dict(
    "TIMESTAMP" => "",
    "RECORD"    => "",
    "WD1"       => "Smp",
    "WD2"       => "Smp",
    "TA"        => "Avg",
    "RH"        => "Avg",
    "HS_Cor"    => "Smp",
    "SFTempK"   => "Avg",
    "SWdown1"   => "Avg",
    "SWdown2"   => "Avg",
    "SWup1"     => "Avg",
    "SWup2"     => "Avg",
    "LWdown1"   => "Avg",
    "LWdown2"   => "Avg",
    "LWup1"     => "Avg",
    "LWup2"     => "Avg",
    "SWdn"      => "Avg",
    "SensorT"   => "Avg",
    "PF_FC4"    => "Avg",
    "WS_FC4"    => "Avg",
    "WS1_Avg"   => "Avg",
    "WS2_Avg"   => "Avg",
    "WS1_Max"   => "Max",
    "WS2_Max"   => "Max",
    "WS1_Std"   => "Std",
    "WS2_Std"   => "Std",
    "Temp_5m_Avg"  => "Avg",
    "RH_5m_Avg"    => "Avg",
    "Temp_16m_Avg" => "Avg",
    "RH_16m_Avg"   => "Avg",
    "Temp_26m_Avg" => "Avg",
    "RH_26m_Avg"   => "Avg",
    "SBTempK"      => "Avg",
    "Ux_Avg"       => "Avg",
    "Uy_Avg"       => "Avg",
    "Uz_Avg"       => "Avg",
    "Ts_Avg"       => "Avg",
    "LI_Pres_Avg"  => "Avg",
    "VP_16m_Avg"   => "Avg",
    "AP_16m_Avg"   => "Avg",
    "Temp_10m_Avg" => "Avg",
    "RH_10m_Avg"   => "Avg",
    "WS_16m_Avg"   => "Avg",
    "WS_16m_Max"   => "Max",
    "WS_16m_Std"   => "Std",
    "WD_16m"       => "Smp",
    "WS_10m_Avg"   => "Avg",
    "WS_10m_Max"   => "Max",
    "WS_10m_Std"   => "Std",
    "WD_10m"       => "Smp",
    "FluxMin_11m"  => "Min",
    "FluxMean_11m" => "Avg",
    "FluxMax_11m"  => "Max",
    "FluxStd_11m"  => "Std",
    "FluxSum_11m"  => "Tot",
    "WindMin_11m"  => "Min",
    "WindMean_11m" => "Avg",
    "WindMax_11m"  => "Max",
    "FC_Flux_11m"  => "Smp",
    "FC_Wind_11m"  => "Smp",
    "FluxMin_4m"   => "Min",
    "FluxMean_4m"  => "Avg",
    "FluxMax_4m"   => "Max",
    "FluxStd_4m"   => "Std",
    "FluxSum_4m"   => "Tot",
    "WindMin_4m"   => "Min",
    "WindMean_4m"  => "Avg",
    "WindMax_4m"   => "Max",
    "FC_Flux_4m"   => "Smp",
    "FC_Wind_4m"   => "Smp",
    "Incoming_SW_16m_Avg"  => "Avg",
    "Outgoing_SW_16m_Avg"  => "Avg",
    "Incoming_LW_16m_Avg"  => "Avg",
    "Outgoing_LW_16m_Avg"  => "Avg",
    "Ux_16m_Avg"   => "Avg",
    "Uy_16m_Avg"   => "Avg",
    "Uz_16m_Avg"   => "Avg",
    "Ts_16m_Avg"   => "Avg",
    "LI_H2Om_16m_Avg" => "Avg",
    "LI_Pres_16m_Avg" => "Avg",
    "VP_26m_Avg"   => "Avg",
    "AP_26m_Avg"   => "Avg",
    "WS_26m_Avg"   => "Avg",
    "WS_26m_Max"   => "Max",
    "WS_26m_Std"   => "Std",
    "WD_26m"       => "Smp",
    "Incoming_SW_26m_Avg"  => "Avg",
    "Outgoing_SW_26m_Avg"  => "Avg",
    "Incoming_UW_26m_Avg"  => "Avg",
    "Outgoing_UW_26m_Avg"  => "Avg",
    "Ux_26m_Avg"   => "Avg",
    "Uy_26m_Avg"   => "Avg",
    "Uz_26m_Avg"   => "Avg",
    "Ts_26m_Avg"   => "Avg",
    "LI_H2Om_26m_Avg" => "Avg",
    "LI_Pres_26m_Avg" => "Avg",
    "WS_5m_Avg"    => "Avg",
    "WS_5m_Max"    => "Max",
    "WS_5m_Std"    => "Std",
    "Ux_5m_Avg"    => "Avg",
    "Uy_5m_Avg"    => "Avg",
    "Uz_5m_Avg"    => "Avg",
    "Ts_5m_Avg"    => "Avg",
    "SWdn_5m"      => "Avg",
    "CS320_T_5m"   => "Avg",
    "wnd_spd_5m"   => "Avg",
)

"""
    save_slow_data(data::DataFrame, output_base::String, sensor::String;
                   nodata::Float64=-9999.0)

Save cleaned slow-frequency data as daily TOA5-style `.dat` files, organized
into monthly subfolders.

Creates files with the naming convention:
    `{output_base}/{YYYYmm}/{sensor}_OneMin_{yyyy-mm-dd}.dat`

Each file has the standard TOA5 4-line header (metadata, column names, units,
processing type) followed by comma-separated data rows.

# Arguments
- `data::DataFrame`: Cleaned slow data with a `TIMESTAMP` column
- `output_base::String`: Root output directory
- `sensor::String`: Sensor name (e.g. "SFC"), used in filenames
- `nodata::Float64`: Value to write for missing data (default: -9999.0)

# Returns
- `Vector{String}`: List of file paths written

# Example
```julia
files = save_slow_data(slow_data, "/path/to/output", "SFC")
```
"""
function save_slow_data(data::DataFrame, output_base::String, sensor::String;
                        nodata::Float64=-9999.0)
    @assert "TIMESTAMP" in names(data) "DataFrame must have a TIMESTAMP column"

    # Get data columns (everything except TIMESTAMP)
    col_names = names(data)
    data_cols = filter(c -> c != "TIMESTAMP", col_names)

    # Group by day
    data.date_ = Date.(data.TIMESTAMP)
    days = sort(unique(data.date_))

    written_files = String[]

    for day in days
        day_data = filter(r -> r.date_ == day, data)

        # Create monthly subfolder: YYYYmm
        month_dir = joinpath(output_base, Dates.format(day, "yyyymm"))
        mkpath(month_dir)

        # Filename: SFC_OneMin_2025-07-15.dat
        fname = joinpath(month_dir, "$(sensor)_OneMin_$(Dates.format(day, "yyyy-mm-dd")).dat")

        open(fname, "w") do io
            # Row 1: TOA5 metadata
            println(io, "\"TOA5\",\"Peddy.jl\",\"processed\",\"$(sensor)\",\"\",\"\",\"\",\"$(sensor)_OneMin\"")

            # Row 2: Column names
            header_cols = ["TIMESTAMP"; data_cols]
            println(io, join(["\"$c\"" for c in header_cols], ","))

            # Row 3: Units
            units = [get(SLOW_DATA_UNITS, c, "") for c in header_cols]
            println(io, join(["\"$u\"" for u in units], ","))

            # Row 4: Processing type
            procs = [get(SLOW_DATA_PROCESSING, c, "Smp") for c in header_cols]
            println(io, join(["\"$p\"" for p in procs], ","))

            # Data rows
            for row in eachrow(day_data)
                vals = String[]
                # Timestamp column
                push!(vals, "\"$(Dates.format(row.TIMESTAMP, "yyyy-mm-dd HH:MM:SS"))\"")
                # Data columns
                for c in data_cols
                    v = row[c]
                    if ismissing(v) || (v isa AbstractFloat && isnan(v))
                        push!(vals, string(nodata))
                    else
                        push!(vals, string(v))
                    end
                end
                println(io, join(vals, ","))
            end
        end

        push!(written_files, fname)
    end

    # Clean up temp column
    select!(data, Not(:date_))

    println("Saved $(length(written_files)) daily files to $output_base")
    return written_files
end


# =============================================================================
# EddyPro external biomet files
# =============================================================================

"""
Column map for EddyPro biomet output.

sensor => ordered vector of `(biomet_name, biomet_unit, source_column, transform)`.
`source_column` may be a vector of candidate names; the first one present in the
DataFrame with at least one non-missing value is used.
Units are spelled the way EddyPro parses them ("C", "%", "kPa", "W+1M-2").
Positional qualifier `_1_1_1` is required by EddyPro's native-header parser.

Notes
- SFC `LI_Pres_Avg` is logged as kPa - 80 (README_output.txt), hence the +80.
- Tower RH columns are already in percent after `clean_tower_slowdata`.
- BOTTOM has no pressure sensor and its RH probe was found stuck (2025), so only
  TA and SW_IN are written; EddyPro falls back to sonic-derived RH and the
  altitude-based pressure for that station.
- In the .eddypro project the biomet columns are referenced by 1-based index
  counting the two timestamp columns: TA=3, RH=4, PA=5, SW_IN=6, LW_IN=7
  (for BOTTOM: TA=3, SW_IN=4).
"""
const BIOMET_COLUMN_MAP = Dict(
    "SFC" => [
        ("TA_1_1_1",    "C",      "TA",          identity),
        ("RH_1_1_1",    "%",      "RH",          identity),
        ("PA_1_1_1",    "kPa",    "LI_Pres_Avg", x -> x + 80.0),
        ("SW_IN_1_1_1", "W+1M-2", "SWdown1",     identity),
        ("LW_IN_1_1_1", "W+1M-2", "LWdown1",     identity),
    ],
    "LOWER" => [
        ("TA_1_1_1",    "C",      "Temp_16m_Avg",        identity),
        ("RH_1_1_1",    "%",      "RH_16m_Avg",          identity),
        ("PA_1_1_1",    "kPa",    "AP_16m_Avg",          identity),
        ("SW_IN_1_1_1", "W+1M-2", "Incoming_SW_16m_Avg", identity),
        ("LW_IN_1_1_1", "W+1M-2", "Incoming_LW_16m_Avg", identity),
    ],
    "UPPER" => [
        ("TA_1_1_1",    "C",      "Temp_26m_Avg",        identity),
        ("RH_1_1_1",    "%",      "RH_26m_Avg",          identity),
        ("PA_1_1_1",    "kPa",    "AP_26m_Avg",          identity),
        ("SW_IN_1_1_1", "W+1M-2", "Incoming_SW_26m_Avg", identity),
        ("LW_IN_1_1_1", "W+1M-2", ["Incoming_LW_26m_Avg", "Incoming_UW_26m_Avg"], identity),  # 2025 logger program labelled it UW
    ],
    "BOTTOM" => [
        ("TA_1_1_1",    "C",      "Temp_5m_Avg", identity),
        ("SW_IN_1_1_1", "W+1M-2", "SWdn_5m",     identity),
    ],
)

# Plausibility bounds per biomet variable; values outside are written as nodata.
const BIOMET_BOUNDS = Dict(
    "TA_1_1_1"    => (-60.0, 20.0),
    "RH_1_1_1"    => (0.0, 100.0),
    "PA_1_1_1"    => (75.0, 95.0),
    "SW_IN_1_1_1" => (-20.0, 1300.0),
    "LW_IN_1_1_1" => (10.0, 400.0),
)

"""
    save_biomet_data(data::DataFrame, output_base::String, sensor::String;
                     nodata::Float64=-9999.0)

Write cleaned slow data as EddyPro external biomet files, one CSV per month:
    `{output_base}/{sensor}_biomet_{yyyymm}.csv`

File layout (what EddyPro 7 expects with "use native header"):
    TIMESTAMP_1,TIMESTAMP_2,TA_1_1_1,RH_1_1_1,PA_1_1_1,SW_IN_1_1_1,LW_IN_1_1_1
    yyyy-mm-dd,HHMM,C,%,kPa,W+1M-2,W+1M-2
    2026-01-02,0000,-5.1200,71.3000,83.0250,412.1000,251.9000

Data are kept at their native (1-min) resolution; EddyPro averages them over
the flux averaging period. Missing or out-of-bounds values become `nodata`.
Column selection and transforms come from `BIOMET_COLUMN_MAP[sensor]`.

# Returns
- `Vector{String}`: list of file paths written
"""
function save_biomet_data(data::DataFrame, output_base::String, sensor::String;
                          nodata::Float64=-9999.0)
    @assert "TIMESTAMP" in names(data) "DataFrame must have a TIMESTAMP column"
    haskey(BIOMET_COLUMN_MAP, sensor) ||
        throw(ArgumentError("No biomet column map defined for sensor \"$sensor\""))
    # Resolve candidate source columns to the one actually present with data
    cmap = [(bname, unit, _resolve_biomet_source(data, src), f)
            for (bname, unit, src, f) in BIOMET_COLUMN_MAP[sensor]]

    missing_src = [bname for (bname, _, src, _) in cmap if src === nothing]
    isempty(missing_src) ||
        @warn "save_biomet_data($sensor): no source column with data for $missing_src, written as nodata"

    mkpath(output_base)
    df = sort(unique(data, :TIMESTAMP), :TIMESTAMP)
    df.month_ = Dates.format.(df.TIMESTAMP, "yyyymm")

    header = join(["TIMESTAMP_1"; "TIMESTAMP_2"; [c[1] for c in cmap]], ",")
    units  = join(["yyyy-mm-dd"; "HHMM"; [c[2] for c in cmap]], ",")

    written = String[]
    for month in sort(unique(df.month_))
        mdf = df[df.month_ .== month, :]
        fname = joinpath(output_base, "$(sensor)_biomet_$(month).csv")
        open(fname, "w") do io
            println(io, header)
            println(io, units)
            for row in eachrow(mdf)
                vals = String[Dates.format(row.TIMESTAMP, "yyyy-mm-dd"),
                              Dates.format(row.TIMESTAMP, "HHMM")]
                for (bname, _, src, f) in cmap
                    v = src === nothing ? missing : row[src]
                    push!(vals, _biomet_value(v, f, get(BIOMET_BOUNDS, bname, (-Inf, Inf)), nodata))
                end
                println(io, join(vals, ","))
            end
        end
        push!(written, fname)
    end

    println("Saved $(length(written)) monthly biomet files to $output_base")
    return written
end

"""Return the first candidate column name that exists in `data` and has any
non-missing, non-NaN value; `nothing` if none does."""
function _resolve_biomet_source(data::DataFrame, src::Union{AbstractString,AbstractVector})
    for c in (src isa AbstractString ? [src] : src)
        c in names(data) || continue
        any(v -> !(ismissing(v) || (v isa AbstractFloat && isnan(v))), data[!, c]) && return String(c)
    end
    return nothing
end

"""Format one biomet value: apply transform, bounds, nodata; 4 decimals."""
function _biomet_value(v, f, bounds, nodata::Float64)
    (ismissing(v) || (v isa AbstractFloat && isnan(v))) && return string(nodata)
    x = f(Float64(v))
    (isnan(x) || x < bounds[1] || x > bounds[2]) && return string(nodata)
    return string(round(x; digits=4))
end
