#!/usr/bin/env julia --project=@.
# Run the Peddy fast-data pipeline for one sensor and one year.
#
# Usage: julia --project=. scripts/run_sensor.jl SENSOR YEAR [--restart] [--out=DIR] [--dates=yyyy-mm-dd,...]
#   SENSOR  SFC | LOWER | UPPER | BOTTOM
#   YEAR    2025 (input 2025/converted) or 2026 (input 2026/data_transfer)
#   --restart  reprocess all dates (existing hourly files are overwritten)
#   --out      processed_HF output folder (default <ROOT>/<YEAR>/processed_HF), e.g. for tests
#   --dates    only these dates
#
# Only dates within YEAR are processed (2025/converted also holds January 2026).
# Slow data are only used as the H2O reference; processed_slow and biomet files
# are not rewritten (2025 biomet files come from eddypro_scripts/make_biomet_2025.py).
# LOWER 2025 has no slow data before 2025-12-27, so it uses the UPPER probe
# (TA and RH together, i.e. the absolute humidity at 25.4 m) for the whole year.

using Peddy
using DimensionalData
using Dates
using DataFrames
using Statistics
using CSV

include("../src/process_sensor.jl")
include("../src/read_data,jl")
include("../src/save_data.jl")

const ROOT = "/home/engbers/Documents/PhD/EC_data_convert"
const INPUT = Dict(2025 => "$ROOT/2025/converted", 2026 => "$ROOT/2026/data_transfer")

opt(name) = (i = findfirst(a -> startswith(a, "--$name="), ARGS); i === nothing ? nothing : String(split(ARGS[i], "=", limit=2)[2]))

sensor = ARGS[1]
year   = parse(Int, ARGS[2])
haskey(SENSORS, sensor) || error("Unknown sensor $sensor")
input_base       = INPUT[year]
processed_output = something(opt("out"), "$ROOT/$year/processed_HF")

# === Slow data (H2O reference) ===
function tower_slow(src)
    suffix = src == "UPPER" ? "_26m" : "_16m"
    raw = read_slow_data(base_path=input_base, sensor=src)
    unique!(raw, :TIMESTAMP)
    sort!(raw, :TIMESTAMP)
    return clean_tower_slowdata(raw, suffix)
end

slow_data = if sensor == "SFC"
    s = clean_slowdata(read_slow_data(base_path=input_base, sensor="SFC"))
    unique!(s, :TIMESTAMP)
    sort!(s, :TIMESTAMP)
    s
elseif sensor == "BOTTOM"
    DataFrame()  # no gas analyzer
elseif sensor == "LOWER" && year == 2025
    println("LOWER 2025: using the UPPER probe as H2O reference")
    up = tower_slow("UPPER")
    select(up, :TIMESTAMP, :Temp_26m_Avg => :Temp_16m_Avg, :RH_26m_Avg => :RH_16m_Avg)
else
    tower_slow(sensor)
end
println("  $(nrow(slow_data)) slow records")

# === Dates ===
dates = if opt("dates") !== nothing
    Date.(split(opt("dates"), ","))
else
    filter(d -> Dates.year(d) == year, list_available_dates(base_path=input_base, sensor=sensor))
end

process_sensor(
    input_base       = input_base,
    sensor           = sensor,
    processed_output = processed_output,
    year             = year,
    dates            = dates,
    resume           = !("--restart" in ARGS),
    slow_data        = slow_data,
)
