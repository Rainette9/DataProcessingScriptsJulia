# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Julia driver around the **Peddy.jl** library for the eddy-covariance (EC) fast data of the Princess
Elisabeth Antarctica (PEA) mast: reads raw TOA5 logger files, runs Peddy's pipeline day by day and writes
hourly processed files plus cleaned 1-min slow data. It is used as a git submodule of `PEA_pipeline`
(github.com/Rainette9/PEA_pipeline), whose cron scripts call `scripts/run_sensor.jl`.

## Architecture

```
scripts/run_sensor.jl            ← The entry point: run_sensor.jl SENSOR [YEAR] [options] (see its header)
scripts/run_mrd.jl               ← MRD post-processor over processed fast data (manual)
scripts/plot_raw_vs_processed.jl ← Raw vs processed plots per day (manual)
src/process_sensor.jl            ← Core: process_sensor(), SENSORS config, make_pipeline(), hourly writer
src/read_data,jl                 ← Readers read_fast_data()/read_slow_data(), slow-data cleaning, H2O reference chain
src/save_data.jl                 ← TOA5-format writer save_slow_data() (OneMin files)
```

**Data flow:** raw TOA5 `.dat` → `read_fast_data()`/`read_slow_data()` → DataFrame → DimArray →
`EddyPipeline.process!()` → `{S}_Fastdata_proc_{date}_{HHMM}.dat` (hourly, ICSV) + `{S}_OneMin_{date}.dat`.
EddyPro biomet files are **not** written here; `PEA_pipeline/make_biomet.py` makes them from the OneMin files.

**Note:** `src/read_data,jl` has a comma in its name, not a period: `include("../src/read_data,jl")`.

## Key Commands

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'                 # once per machine
julia --project=. scripts/run_sensor.jl LOWER                        # resume mode (daily cron)
julia --project=. scripts/run_sensor.jl LOWER --restart --from 2026-10-20 --to 2026-10-20
```

Paths come from the `PEA_*` environment (`PEA_INPUT_BASE`, `PEA_HF_OUT`, `PEA_LF_OUT`, `PEA_YEAR`) or the
`--input= --out= --slow_out=` options; the `EC_ROOT` layout is the PC alternative. All documented in the
header of `scripts/run_sensor.jl`.

Julia 1.11+ is required. Peddy.jl is pinned in `Manifest.toml` to branch `fix_licor` of
github.com/pleibers/Peddy.jl at the tree of commit fe0c4c3 (the "v3" processing: sonic `diag_margin`,
Ts rolling-median check, `couple_sonic` despiking). `Pkg.update("Peddy")` would move it; do not run it
without re-validating the output. There is no test suite here; validate changes by reprocessing one day into
a scratch folder and comparing with `PEA_pipeline/tools/compare_hf.py` (0 differing values expected unless
the change is meant to alter the output).

## Processing versions

- **v3** (reference since Oct 2026, identical on all machines): `SENSORS` in `src/process_sensor.jl`
  (`diag_margin` SFC 0 / towers 2, `pres_offset` 80 kPa for SFC), `ts_max_deviation` 10 K,
  `couple_sonic=true`, `max_gap_samples` 5, H2O bounds 0–400 mmol/m³, reference grid with ≤10 min
  interpolation. Overrides: `--diag_margin=`, `--ts_max_deviation=` on `run_sensor.jl`.
- The year only selects the LI-7500 calibration coefficients (`Peddy.default_calibration_coefficients`,
  SFC 2024 vs ≥2025) and filters the dates; 2025 (logger-card raw, `TOA5_STN1fast_<N>.dat`) and 2026+
  (FTP hourly files merged to daily) go through the same code. 2025 quirks are detected from the data
  (RH fraction → percent, missing own slow data → reference chain).

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

`clean_tower_slowdata` returns RH in **percent** (raw tower RH_16m/RH_26m are 0–1 and are auto-scaled; RH_10m/RH_5m are already percent); `build_lf_dimarray`/`_lf_grid` divide by 100 before handing the reference to Peddy. Logger code `-7999` is masked there as well.

## H2O reference chain (LOWER/UPPER)

`H2O_REF_CHAIN` in `src/read_data,jl`: LOWER = own 16 m → UPPER 26 m → LOWER-logger 10 m probe; UPPER = own 26 m → 10 m probe. `build_reference_chain` merges the sources minute by minute (first valid wins, `REF_SRC` = source index, 9 = interpolated gap ≤ 10 min). `process_sensor(...; ref_data=ref)` uses it instead of `slow_data`; the reference is placed on a regular 1-min grid from day−1 h to day+1 d+1 h with NaN gaps, so Peddy's nearest-minute matching never uses a reference hours away (NaN reference → NaN H2O). After the pipeline H2O outside 0–400 mmol/m³ is set to NaN. Per-day source counts go to `{S}_h2o_ref_{YYYYMM}.csv` next to the processing log. Only the sensor's own cleaned slow data is written as OneMin files. Days without any reference minute get NaN H2O. SFC uses its own cleaned slow data (no chain); BOTTOM has no gas analyzer.

## Peddy.jl Pipeline Steps

Configured via `EddyPipeline(; sensor, quality_control, gas_analyzer, despiking, make_continuous, gap_filling, double_rotation, mrd, output)`. Each step is optional (set to `nothing` to skip). Execution order is fixed: QC → gas analyzer → despiking → make continuous → gap filling → double rotation → MRD → output.

## DimArray Conventions

- High-frequency data: dimensions `(Var, Ti)` — variables as rows, time as columns
- **ICSVOutput requires `(Ti, Var)` orientation** — transpose before writing with `hcat()` + `DimArray(matrix, (Ti(times), Var(vars)))`

## Sensor Data

Four stations: `"SFC"`, `"LOWER"`, `"UPPER"`, `"BOTTOM"`. TOA5 format: 4 header rows (metadata, column names, units, processing type), then CSV data. The readers match files by name (`fast`, `OneMin`/`slow`, `wind`) inside any folder whose path contains the sensor name, and take dates from the TIMESTAMP column, so both the card and the FTP naming work.

## Known Pitfalls

- CSAT3 diagnostic (`diag_csat`) is a 6-bit counter (0–63) + error flags in bits 6+. Preprocess with `Float64.(fast_data.diag_csat .>= 64)` so only actual errors are flagged.
- LI-7500 diagnostic: if constant at 247 (above default threshold 240), set `diag_gas=255` in `LICOR()` constructor to avoid NaN-ing all H2O/P data.
- `skipmissing()` does NOT filter `Float64` NaN values — use `filter(!isnan, ...)` for statistics on pipeline output.
- ICSV output uses `nodata=-9999.0`; read back with `CSV.read(...; missingstring="-9999.0")`.
- Raw daily fast files span ~25 h and `read_fast_data(dates=[d])` also loads the neighbouring file; `process_sensor` filters each day to its calendar day before processing (otherwise ~2 days are processed and the neighbour's hourly files overwritten).
- `process_sensor` appends to an existing monthly processing log when `resume=true` (the daily cron) and only starts a fresh log with `--restart`.
- `data/_raw/{SENSOR}` is walked by the readers too (path contains the sensor name); harmless but slow.
