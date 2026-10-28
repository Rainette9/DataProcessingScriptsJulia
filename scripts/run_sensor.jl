#!/usr/bin/env julia
# Parameterized EC processing driver for one sensor (used by run_LOWER.jl / run_UPPER.jl).
#
# Usage:
#   julia --project=<DataProcessingScriptsJulia> scripts/run_sensor.jl SENSOR [options]
#     --restart        reprocess all dates (default: resume, skip dates with output)
#     --from DATE      first date to process (YYYY-MM-DD)
#     --to DATE        last date to process (YYYY-MM-DD)
#     --slow-only      only read/clean/save the slow data, no fast-data processing
#     --no-save-slow   do not (re)write the OneMin slow-data files
#
# Environment (defaults = production paths used by the daily cron):
#   PEA_INPUT_BASE   raw data root        (default /data/disk1/shared/PEA/data)
#   PEA_HF_OUT       processed_HF root    (default /data/disk1/shared/PEA/processed_HF)
#   PEA_LF_OUT       processed_LF root    (default /data/disk1/shared/PEA/processed_LF; sensor appended)
#   PEA_YEAR         year (calibration coefficients + only dates of this year are processed)
#   PEA_REF_CHAIN_<SENSOR>  override of the H2O reference chain, "SENSOR:suffix:TAcol:RHcol,..."
#   PEA_GAP_FILLING  0 = no interpolation of short gaps (flagged/despiked samples stay NaN; default 1)
#
# H2O reference: see H2O_REF_CHAIN in src/read_data,jl (own probe → UPPER 26 m → LOWER 10 m
# for LOWER; own probe → LOWER 10 m for UPPER). Only the sensor's own cleaned slow data is
# saved as OneMin files; the reference itself is never written as slow data.

using Peddy
using DimensionalData
using Dates
using DataFrames
using Statistics
using CSV

include("../src/process_sensor.jl")
include("../src/read_data,jl")
include("../src/save_data.jl")

# === Arguments ===
isempty(ARGS) && error("usage: run_sensor.jl SENSOR [--restart] [--from YYYY-MM-DD] [--to YYYY-MM-DD] [--slow-only] [--no-save-slow]")
sensor = String(ARGS[1])
sensor in ("LOWER", "UPPER", "BOTTOM") || error("run_sensor.jl supports tower sensors (LOWER/UPPER/BOTTOM), got $sensor")

function argval(flag)
    i = findfirst(==(flag), ARGS)
    i === nothing && return nothing
    i < length(ARGS) || error("$flag needs a value")
    return ARGS[i + 1]
end
restart      = "--restart" in ARGS
slow_only    = "--slow-only" in ARGS
save_slow    = !("--no-save-slow" in ARGS)
date_from    = (v = argval("--from")) === nothing ? nothing : Date(v)
date_to      = (v = argval("--to"))   === nothing ? nothing : Date(v)

# === Configuration ===
input_base       = get(ENV, "PEA_INPUT_BASE", "/data/disk1/shared/PEA/data")
processed_output = get(ENV, "PEA_HF_OUT",     "/data/disk1/shared/PEA/processed_HF")
slow_output      = joinpath(get(ENV, "PEA_LF_OUT", "/data/disk1/shared/PEA/processed_LF"), sensor)
year             = parse(Int, get(ENV, "PEA_YEAR", "2026"))
own_suffix       = SENSORS[sensor].suffix
gap_filling      = get(ENV, "PEA_GAP_FILLING", "1") != "0"

chain_spec = get(ENV, "PEA_REF_CHAIN_$(sensor)", "")
chain = isempty(chain_spec) ? H2O_REF_CHAIN[sensor] : parse_ref_chain(chain_spec)

println("=" ^ 60)
println("run_sensor.jl  sensor=$sensor  year=$year")
println("  input_base       = $input_base")
println("  processed_output = $processed_output")
println("  slow_output      = $slow_output")
println("  restart=$restart  slow_only=$slow_only  save_slow=$save_slow  from=$date_from  to=$date_to  gap_filling=$gap_filling")
println("  reference chain  = ", join(["$(e.sensor):$(e.ta)/$(e.rh)" for e in chain], " → "))
println("=" ^ 60)

# === Slow data: own sensor (saved) + reference chain (not saved) ===
cache = Dict{String, DataFrame}()
slow_data = load_cleaned_slow!(cache, sensor, own_suffix; input_base=input_base)

if save_slow && nrow(slow_data) > 0
    own_year = filter(r -> Dates.year(r.TIMESTAMP) == year, slow_data)
    println("Saving slow data (OneMin) for $sensor, year $year: $(nrow(own_year)) records")
    nrow(own_year) > 0 && save_slow_data(own_year, slow_output, sensor)
end

ref_data = build_reference_chain(chain; input_base=input_base, cache=cache)
nrow(ref_data) > 0 && interpolate_short_gaps!(ref_data; max_gap_minutes=10)

if slow_only
    println("--slow-only: done.")
    exit(0)
end

# === Dates ===
dates = list_available_dates(base_path=input_base, sensor=sensor)
dates = filter(d -> Dates.year(d) == year, dates)
date_from === nothing || (dates = filter(d -> d >= date_from, dates))
date_to   === nothing || (dates = filter(d -> d <= date_to, dates))
println("Dates selected: $(length(dates))" * (isempty(dates) ? "" : " ($(first(dates)) to $(last(dates)))"))
if isempty(dates)
    println("Nothing to process.")
    exit(0)
end

# === Process ===
process_sensor(
    input_base       = input_base,
    sensor           = sensor,
    processed_output = processed_output,
    year             = year,
    dates            = dates,
    resume           = !restart,
    slow_data        = slow_data,
    ref_data         = nrow(ref_data) > 0 ? ref_data : nothing,
    config           = Dict(:gap_filling => gap_filling),
)
