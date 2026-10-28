#!/usr/bin/env julia --project=@.
# Read, clean and save the 1-min slow data of one sensor and write the monthly
# EddyPro biomet files (TA, RH, PA, SW_IN, LW_IN; see BIOMET_COLUMN_MAP in src/save_data.jl).
#
# Usage: julia --project=. scripts/run_slow.jl SENSOR YEAR [--input=DIR] [--slow_out=DIR] [--biomet_out=DIR]
#   defaults: input as in run_sensor.jl, slow_out <ROOT>/<YEAR>/processed_slow/<SENSOR>,
#             biomet_out <ROOT>/<YEAR>/biomet/<SENSOR>; ROOT from the environment variable EC_ROOT.
# Only records within YEAR are written.

using Dates
using DataFrames
using CSV
using Statistics

include("../src/read_data,jl")
include("../src/save_data.jl")

const ROOT = get(ENV, "EC_ROOT", "/home/engbers/Documents/PhD/EC_data_convert")
const INPUT = Dict(2025 => "$ROOT/2025/converted", 2026 => "$ROOT/2026/data_transfer")
const SUFFIX = Dict("LOWER" => "_16m", "UPPER" => "_26m", "BOTTOM" => "_5m")

opt(name) = (i = findfirst(a -> startswith(a, "--$name="), ARGS); i === nothing ? nothing : String(split(ARGS[i], "=", limit=2)[2]))

sensor = ARGS[1]
year   = parse(Int, ARGS[2])
input_base = something(opt("input"), get(INPUT, year, nothing), "")
isdir(input_base) || error("No raw data folder for $year: pass --input=DIR")
slow_out   = something(opt("slow_out"), "$ROOT/$year/processed_slow/$sensor")
biomet_out = something(opt("biomet_out"), "$ROOT/$year/biomet/$sensor")

println("Reading slow data for $sensor from $input_base")
raw = read_slow_data(base_path=input_base, sensor=sensor)
slow = if sensor == "SFC"
    clean_slowdata(raw)
else
    unique!(raw, :TIMESTAMP)
    sort!(raw, :TIMESTAMP)
    clean_tower_slowdata(raw, SUFFIX[sensor])  # RH fraction -> percent
end
unique!(slow, :TIMESTAMP)
sort!(slow, :TIMESTAMP)
filter!(r -> Dates.year(r.TIMESTAMP) == year, slow)
println("  $(nrow(slow)) slow records in $year")

save_slow_data(slow, slow_out, sensor)
save_biomet_data(slow, biomet_out, sensor)
