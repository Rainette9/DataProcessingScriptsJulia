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
