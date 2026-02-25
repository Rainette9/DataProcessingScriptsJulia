#!/usr/bin/env julia --project=.
# Run the EC processing pipeline for the SFC sensor.

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
input_base       = "/home/engbers/Documents/PhD/EC_data_convert/2025/converted"
processed_output = "/home/engbers/Documents/PhD/EC_data_convert/2025/processed_HF"
year             = 2025

# === Load and clean slow data ===
println("Reading slow data for SFC...")
slow_raw  = read_slow_data(base_path=input_base, sensor="SFC")
slow_data = clean_slowdata(slow_raw)
unique!(slow_data, :TIMESTAMP)
sort!(slow_data, :TIMESTAMP)
println("  $(nrow(slow_data)) slow records loaded")

# === Save slow data (including wind columns) ===
slow_output = "/home/engbers/Documents/PhD/EC_data_convert/2025/processed_slow/SFC"
println("Saving slow data (OneMin) with wind columns...")
save_slow_data(slow_data, slow_output, "SFC")

# === Process ===
# Pass --restart to reprocess all dates; default resumes from where it left off.
process_sensor(
    input_base       = input_base,
    sensor           = "SFC",
    processed_output = processed_output,
    year             = year,
    resume           = !("--restart" in ARGS),
    slow_data        = slow_data,
)
