#!/usr/bin/env julia --project=.
# Run the EC processing pipeline for the LOWER sensor.

using Peddy
using DimensionalData
using Dates
using DataFrames
using Statistics
using CSV

include("process_sensor.jl")
include("src/read_data,jl")
include("save_data.jl")

# === Configuration ===
input_base       = "/home/engbers/Documents/PhD/EC_data_convert/2026/data_transfer"
processed_output = "/home/engbers/Documents/PhD/EC_data_convert/2026/processed_HF"
year             = 2026

# === Load slow data ===
println("Reading slow data for LOWER...")
slow_raw = read_slow_data(base_path=input_base, sensor="LOWER")
unique!(slow_raw, :TIMESTAMP)
sort!(slow_raw, :TIMESTAMP)
slow_data = clean_tower_slowdata(slow_raw, "_16m")
println("  $(nrow(slow_data)) slow records loaded")

# === Save slow data ===
slow_output = "/home/engbers/Documents/PhD/EC_data_convert/2025/processed_slow/LOWER"
println("Saving slow data (OneMin)...")
save_slow_data(slow_data, slow_output, "LOWER")

# === Process ===
# Pass --restart to reprocess all dates; default resumes from where it left off.
process_sensor(
    input_base       = input_base,
    sensor           = "LOWER",
    processed_output = processed_output,
    year             = year,
    resume           = !("--restart" in ARGS),
    slow_data        = slow_data,
)
