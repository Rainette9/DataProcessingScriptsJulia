# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Eddy covariance (EC) data processing pipeline for snow research stations, built on the **Peddy.jl** library. The main workflow runs in `EC_processing.ipynb` (Jupyter notebook), reading raw TOA5 sensor data, processing it through Peddy's pipeline, and outputting results.

## Architecture

```
EC_processing.ipynb              ← Main pipeline notebook (runs end-to-end)
src/read_data,jl                 ← Data readers: read_fast_data(), read_slow_data(), clean_slowdata()
src/save_data.jl                 ← TOA5-format writer: save_slow_data()
src/process_sensor.jl            ← Core pipeline engine: process_sensor()
scripts/run_SFC.jl               ← Per-station entry points (SFC/BOTTOM: own slow data)
scripts/run_BOTTOM.jl
scripts/run_sensor.jl            ← Env-configurable driver with the H2O reference chain (LOWER/UPPER)
scripts/run_LOWER.jl             ← Wrappers around run_sensor.jl
scripts/run_UPPER.jl
scripts/run_mrd.jl               ← MRD post-processor
scripts/plot_raw_vs_processed.jl ← Visualization
Peddy.jl/                        ← Local dev copy of the EC processing library (separate repo)
```

**Data flow:** Raw TOA5 `.dat` files → `read_fast_data()`/`read_slow_data()` → DataFrame → convert to DimArray → `EddyPipeline.process!()` → output (ICSV/NetCDF/Memory)

**Note:** The file `src/read_data,jl` has a comma in its name, not a period. Include it with `include("../src/read_data,jl")` from scripts, or `include("src/read_data,jl")` from the project root.

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

`clean_tower_slowdata` returns RH in **percent** (raw tower RH_16m/RH_26m are 0–1 and are auto-scaled; RH_10m/RH_5m are already percent); `build_lf_dimarray`/`_regular_lf_window` divide by 100 before handing the reference to Peddy. Logger code `-7999` is masked there as well.

## H2O reference chain (LOWER/UPPER)

`H2O_REF_CHAIN` in `src/read_data,jl`: LOWER = own 16 m → UPPER 26 m → LOWER-logger 10 m probe; UPPER = own 26 m → 10 m probe. `build_reference_chain` merges the sources minute by minute (first valid wins, `REF_SRC` = source index, 9 = interpolated gap ≤ 10 min). `process_sensor(...; ref_data=ref)` uses it instead of `slow_data`; the reference is placed on a regular 1-min grid spanning the day's fast data with NaN gaps, so Peddy's nearest-minute matching never uses a reference hours away (NaN reference → NaN H2O). After the pipeline H2O outside 0–400 mmol/m³ is set to NaN. Per-day source counts go to `{S}_h2o_ref_{YYYYMM}.csv` next to the processing log. Only the sensor's own cleaned slow data is written as OneMin files.

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
- Raw daily fast files span ~25 h and `read_fast_data(dates=[d])` also loads the neighbouring file; `process_sensor` filters each day to its calendar day before processing (otherwise ~2 days are processed and the neighbour's hourly files overwritten).
- `process_sensor` appends to an existing monthly processing log when `resume=true` (the daily cron) and only starts a fresh log with `--restart`.
- `data/_raw/{SENSOR}` is walked by the readers too (path contains the sensor name); harmless but slow.
