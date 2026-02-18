#!/usr/bin/env julia --project=.
# Run the EC processing pipeline for the BOTTOM sensor (CSAT3 only, no LICOR).

using Peddy
using DimensionalData
using Dates
using DataFrames
using Statistics
using CSV

include("process_sensor.jl")
include("src/read_data,jl")

# === Configuration ===
input_base       = "/home/engbers/Documents/PhD/EC_data_convert/2025/converted"
processed_output = "/home/engbers/Documents/PhD/EC_data_convert/2025/processed_HF"
year             = 2025

# === Load slow data ===
println("Reading slow data for BOTTOM...")
slow_data = read_slow_data(base_path=input_base, sensor="BOTTOM")
unique!(slow_data, :TIMESTAMP)
sort!(slow_data, :TIMESTAMP)
println("  $(nrow(slow_data)) slow records loaded")

# === Process ===
# Pass --restart to reprocess all dates; default resumes from where it left off.
process_sensor(
    input_base       = input_base,
    sensor           = "BOTTOM",
    processed_output = processed_output,
    year             = year,
    resume           = !("--restart" in ARGS),
    slow_data        = slow_data,
)
