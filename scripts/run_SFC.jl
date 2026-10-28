#!/usr/bin/env julia
# Run the EC processing pipeline for the SFC sensor (own slow data as H2O reference).
# Thin wrapper around run_sensor.jl: resume mode by default; --restart / --from / --to /
# --slow-only / --no-save-slow and the PEA_* environment are documented there.
pushfirst!(ARGS, "SFC")
include("run_sensor.jl")
