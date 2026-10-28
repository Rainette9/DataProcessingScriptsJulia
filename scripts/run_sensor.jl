#!/usr/bin/env julia
# Peddy fast-data driver for one sensor (SFC, LOWER, UPPER or BOTTOM).
# The thin wrappers run_SFC.jl / run_LOWER.jl / run_UPPER.jl / run_BOTTOM.jl (daily cron,
# resume mode) include this file; reruns call it directly.
#
# Usage: julia --project=<DataProcessingScriptsJulia> scripts/run_sensor.jl SENSOR [YEAR] [options]
#   --restart                 reprocess all dates (default: resume, skip dates with output)
#   --from DATE, --to DATE    date range (YYYY-MM-DD); --from=DATE also works
#   --dates=d1,d2,...         only these dates
#   --months=yyyymm,...       only these months (days without raw data are skipped)
#   --slow-only               only read/clean/save the slow data (OneMin), no fast-data processing
#   --no-save-slow            do not (re)write the OneMin slow-data files
#   --input=DIR, --out=DIR, --slow_out=DIR   override the paths below
#   --diag_margin=N, --ts_max_deviation=K    override process_sensor defaults (SFC 0 / towers 2; 10 K)
#   --logtag=SUFFIX           suffix for the monthly processing log names (partial runs)
#
# Paths, first match wins: option, PEA_* environment, EC_ROOT layout, workstation default.
#   PEA_INPUT_BASE  raw data root (default /data/disk1/shared/PEA/data;
#                   EC_ROOT: <EC_ROOT>/2025/converted, <EC_ROOT>/<YEAR>/data_transfer otherwise)
#   PEA_HF_OUT      processed_HF root (default /data/disk1/shared/PEA/processed_HF; EC_ROOT: <EC_ROOT>/<YEAR>/processed_HF)
#   PEA_LF_OUT      processed_LF root, sensor appended (default /data/disk1/shared/PEA/processed_LF;
#                   EC_ROOT: <EC_ROOT>/<YEAR>/processed_slow)
#   PEA_YEAR        year (default: YEAR argument, else the current year); selects the calibration
#                   coefficients and restricts the dates to that year
#   PEA_REF_CHAIN_<SENSOR>  H2O reference chain override, "SENSOR:suffix:TAcol:RHcol,..."
#   PEA_GAP_FILLING 0 = no interpolation of short gaps (flagged/despiked samples stay NaN; default 1)
#
# H2O reference: SFC uses its own cleaned slow data (TA/RH); LOWER/UPPER use H2O_REF_CHAIN in
# src/read_data,jl (own probe → UPPER 26 m → LOWER-logger 10 m for LOWER; own probe → 10 m for
# UPPER); BOTTOM has no gas analyzer. Only the sensor's own cleaned slow data is written as OneMin
# files; the reference itself is never written. EddyPro biomet files are made from the OneMin
# files by transfer_scripts/make_biomet.py (save_biomet_data in src/save_data.jl is the
# equivalent Julia writer; it is not called here).

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
const USAGE = "usage: run_sensor.jl SENSOR [YEAR] [--restart] [--from D] [--to D] [--dates=..] [--months=..] " *
              "[--slow-only] [--no-save-slow] [--input=DIR] [--out=DIR] [--slow_out=DIR] " *
              "[--diag_margin=N] [--ts_max_deviation=K] [--logtag=S]"
isempty(ARGS) && error(USAGE)
sensor = String(ARGS[1])
haskey(SENSORS, sensor) || error("Unknown sensor $sensor (SFC, LOWER, UPPER, BOTTOM)\n$USAGE")
year_arg = length(ARGS) >= 2 && occursin(r"^\d{4}$", ARGS[2]) ? parse(Int, ARGS[2]) : nothing

# "--name=value" options
opt(name) = (i = findfirst(a -> startswith(a, "--$name="), ARGS); i === nothing ? nothing : String(split(ARGS[i], "=", limit=2)[2]))
# "--flag value" options (also accept "--flag=value")
function argval(flag)
    v = opt(flag[3:end])
    v === nothing || return v
    i = findfirst(==(flag), ARGS)
    i === nothing && return nothing
    i < length(ARGS) || error("$flag needs a value")
    return ARGS[i + 1]
end
envval(k) = (v = get(ENV, k, ""); isempty(v) ? nothing : v)

restart      = "--restart" in ARGS
slow_only    = "--slow-only" in ARGS
save_slow    = !("--no-save-slow" in ARGS)
date_from    = (v = argval("--from")) === nothing ? nothing : Date(v)
date_to      = (v = argval("--to"))   === nothing ? nothing : Date(v)

# === Configuration ===
year    = something(year_arg, (v = envval("PEA_YEAR"); v === nothing ? nothing : parse(Int, v)), Dates.year(today()))
ec_root = envval("EC_ROOT")
ec_in   = ec_root === nothing ? nothing : (year == 2025 ? "$ec_root/2025/converted" : "$ec_root/$year/data_transfer")
ec_hf   = ec_root === nothing ? nothing : "$ec_root/$year/processed_HF"
ec_lf   = ec_root === nothing ? nothing : "$ec_root/$year/processed_slow"

input_base       = something(opt("input"),    envval("PEA_INPUT_BASE"), ec_in, "/data/disk1/shared/PEA/data")
processed_output = something(opt("out"),      envval("PEA_HF_OUT"),     ec_hf, "/data/disk1/shared/PEA/processed_HF")
slow_output      = joinpath(something(opt("slow_out"), envval("PEA_LF_OUT"), ec_lf, "/data/disk1/shared/PEA/processed_LF"), sensor)
isdir(input_base) || error("Raw data folder not found: $input_base (pass --input=DIR or set PEA_INPUT_BASE)")
own_suffix       = SENSORS[sensor].suffix
gap_filling      = get(ENV, "PEA_GAP_FILLING", "1") != "0"
use_chain        = haskey(H2O_REF_CHAIN, sensor)

chain = if use_chain
    spec = get(ENV, "PEA_REF_CHAIN_$(sensor)", "")
    isempty(spec) ? H2O_REF_CHAIN[sensor] : parse_ref_chain(spec)
else
    nothing
end

config = Dict{Symbol, Any}(:gap_filling => gap_filling)
opt("diag_margin")      === nothing || (config[:diag_margin]      = parse(Int, opt("diag_margin")))
opt("ts_max_deviation") === nothing || (config[:ts_max_deviation] = parse(Float64, opt("ts_max_deviation")))

println("=" ^ 60)
println("run_sensor.jl  sensor=$sensor  year=$year")
println("  input_base       = $input_base")
println("  processed_output = $processed_output")
println("  slow_output      = $slow_output")
println("  restart=$restart  slow_only=$slow_only  save_slow=$save_slow  from=$date_from  to=$date_to  gap_filling=$gap_filling")
use_chain && println("  reference chain  = ", join(["$(e.sensor):$(e.ta)/$(e.rh)" for e in chain], " → "))
isempty(config) || println("  config           = $config")
println("=" ^ 60)

# === Slow data: own sensor (saved as OneMin) + reference (not saved) ===
cache = Dict{String, DataFrame}()
slow_data = if sensor == "SFC"
    raw = read_slow_data(base_path=input_base, sensor="SFC")
    if nrow(raw) == 0
        @warn "No slow data found for SFC under $input_base"
        DataFrame(TIMESTAMP=DateTime[])
    else
        s = clean_slowdata(raw)
        unique!(s, :TIMESTAMP)
        sort!(s, :TIMESTAMP)
        println("  $(nrow(s)) slow records for SFC")
        s
    end
else
    load_cleaned_slow!(cache, sensor, own_suffix; input_base=input_base)
end

if save_slow && nrow(slow_data) > 0
    own_year = filter(r -> Dates.year(r.TIMESTAMP) == year, slow_data)
    println("Saving slow data (OneMin) for $sensor, year $year: $(nrow(own_year)) records")
    nrow(own_year) > 0 && save_slow_data(own_year, slow_output, sensor)
end

ref_data = nothing
if use_chain
    ref = build_reference_chain(chain; input_base=input_base, cache=cache)
    if nrow(ref) > 0
        interpolate_short_gaps!(ref; max_gap_minutes=MAX_REF_GAP_MIN)
        ref_data = ref
    end
end

if slow_only
    println("--slow-only: done.")
    exit(0)
end

# === Dates ===
dates = if opt("dates") !== nothing
    Date.(split(opt("dates"), ","))
elseif opt("months") !== nothing
    [d for m in split(opt("months"), ",")
       for d in Date(parse(Int, m[1:4]), parse(Int, m[5:6]), 1):Day(1):lastdayofmonth(Date(parse(Int, m[1:4]), parse(Int, m[5:6]), 1))]
else
    filter(d -> Dates.year(d) == year, list_available_dates(base_path=input_base, sensor=sensor))
end
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
    ref_data         = ref_data,
    log_tag          = something(opt("logtag"), ""),
    config           = config,
)
