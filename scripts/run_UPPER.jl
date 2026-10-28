#!/usr/bin/env julia
# Run the EC processing pipeline for the UPPER sensor (H2O reference chain).
# Thin wrapper around run_sensor.jl: resume mode by default; --restart / --from / --to /
# --slow-only / --no-save-slow and the PEA_* environment are documented there.
pushfirst!(ARGS, "UPPER")
include("run_sensor.jl")
