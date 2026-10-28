#!/usr/bin/env julia
# Run the EC processing pipeline for the UPPER sensor.
# Thin wrapper around run_sensor.jl (reference chain for the H2O correction, env-configurable
# paths, --restart / --from / --to / --slow-only / --no-save-slow). See scripts/run_sensor.jl.
pushfirst!(ARGS, "UPPER")
include("run_sensor.jl")
