#!/usr/bin/env julia
# Run the EC processing pipeline for the BOTTOM sensor (CSAT3 only, no LICOR).
# Thin wrapper around run_sensor.jl: resume mode by default; --restart / --from / --to /
# --slow-only / --no-save-slow and the PEA_* environment are documented there.
pushfirst!(ARGS, "BOTTOM")
include("run_sensor.jl")
