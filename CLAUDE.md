# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Eddy covariance (EC) data processing pipeline for snow research stations, built on the **Peddy.jl** library. The main workflow runs in `EC_processing.ipynb` (Jupyter notebook), reading raw TOA5 sensor data, processing it through Peddy's pipeline, and outputting results.

## Architecture

```
EC_processing.ipynb          ← Main pipeline notebook (runs end-to-end)
src/read_data,jl             ← Data readers: read_fast_data(), read_slow_data(), clean_slowdata()
save_data.jl                 ← TOA5-format writer: save_slow_data()
Peddy.jl/                    ← Local dev copy of the EC processing library (separate repo)
```

**Data flow:** Raw TOA5 `.dat` files → `read_fast_data()`/`read_slow_data()` → DataFrame → convert to DimArray → `EddyPipeline.process!()` → output (ICSV/NetCDF/Memory)

**Note:** The file `src/read_data,jl` has a comma in its name, not a period. Include it with `include("src/read_data,jl")`.

## Key Commands

```bash
# Activate project environment and run notebook
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run Peddy.jl tests
julia --project=Peddy.jl -e 'using Pkg; Pkg.test()'

# Build Peddy.jl docs
julia --project=Peddy.jl/docs -e 'using Pkg; Pkg.instantiate(); include("Peddy.jl/docs/make.jl")'
```

Julia 1.11+ is required.

## Peddy.jl Unit Conventions (Critical)

These are verified from source and diverge from what some metadata files claim:

| Variable | Unit | Notes |
|----------|------|-------|
| Ux, Uy, Uz | m/s | |
| Ts | °C | Despite `variable_metadata.jl` saying "K", bounds check is -50 to 50 |
| H2O | mmol/m³ | LI-7500 raw output; H2OCalibration corrects it |
| P | kPa | Peddy multiplies by 1000 internally when Pa needed |
| TA (slow) | °C | Magnus formula expects Celsius |
| RH (slow) | fraction 0–1 | Peddy does `* 100` internally — do NOT pass percentage |

## Peddy.jl Pipeline Steps

Configured via `EddyPipeline(; sensor, quality_control, gas_analyzer, despiking, make_continuous, gap_filling, double_rotation, mrd, output)`. Each step is optional (set to `nothing` to skip). Execution order is fixed: QC → gas analyzer → despiking → make continuous → gap filling → double rotation → MRD → output.

## DimArray Conventions

- High-frequency data: dimensions `(Var, Ti)` — variables as rows, time as columns
- **ICSVOutput requires `(Ti, Var)` orientation** — transpose before writing with `hcat()` + `DimArray(matrix, (Ti(times), Var(vars)))`

## Sensor Data

Four stations: `"SFC"`, `"LOWER"`, `"UPPER"`, `"BOTTOM"`. TOA5 format: 4 header rows (metadata, column names, units, processing type), then CSV data.

## Known Pitfalls

- CSAT3 diagnostic (`diag_csat`) is a 6-bit counter (0–63) + error flags in bits 6+. Preprocess with `Float64.(fast_data.diag_csat .>= 64)` so only actual errors are flagged.
- LI-7500 diagnostic: if constant at 247 (above default threshold 240), set `diag_gas=255` in `LICOR()` constructor to avoid NaN-ing all H2O/P data.
- `skipmissing()` does NOT filter `Float64` NaN values — use `filter(!isnan, ...)` for statistics on pipeline output.
- ICSV output uses `nodata=-9999.0`; read back with `CSV.read(...; missingstring="-9999.0")`.
