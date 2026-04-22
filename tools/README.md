# C-GEM Analysis Tools

Python scripts for visualizing, validating, and analyzing C-GEM model output.
Each script produces **one comprehensive, publication-quality figure**.

## Requirements

```bash
pip install -r tools/requirements.txt
```

All scripts must be run **from the project root directory** (where `OUT/` is located).

`check_inputs.py` and `convert_obs_to_model_units.py` can be used before the
model is run. The plotting scripts require model outputs in `OUT/`.

## Design Principles

1. **One figure per script** — comprehensive multi-panel layout, publication-quality (dpi 250).
2. **All tools work without observation data** — they plot model output alone.
3. **Observation overlay is optional** — pass `--with-obs` (uses default files in
   `INPUT/Validation/`) or `--obs-file PATH` for custom data.
4. **Observation data should be supplied in model units**.
   The bundled compact validation tables already follow this convention; see
   `INPUT/Validation/README.md` for the conversion reference.
5. **No hardcoded paths** — all file paths are command-line arguments.
6. **Universal** — not tied to any specific estuary.
7. **Timestamped outputs** — every run saves a dated copy + a "latest" copy.

## Quick Start

After running the model (`bin/Release/C_GEM_Daily`), the output is in `OUT/`.
The helpers read both the default CSV files and CGEM binary (`.bin`) files.

```bash
# ---- 1. Pre-run: check input data quality ----
python tools/check_inputs.py
# → OUT/figures/fig_input_diagnostics.png  +  OUT/tables/check_inputs_report.csv

# ---- 1b. Optional: convert field-unit observation tables to model-unit copies ----
python tools/convert_obs_to_model_units.py
# → INPUT/Validation/converted/*.csv

# ---- 2. Longitudinal WQ profiles (dry/wet bands + obs scatter + metrics) ----
python tools/plot_profiles.py                     # model-only
python tools/plot_profiles.py --with-obs          # + obs overlay + RMSE/R²
python tools/plot_profiles.py --vars Sal O2 Phy1  # specific variables
# → OUT/figures/fig_profiles.png  +  OUT/tables/validation_profiles.csv

# ---- 3. Seasonal time-series at monitoring stations ----
python tools/plot_seasonal.py                     # model-only (auto-detect stations)
python tools/plot_seasonal.py --with-obs          # + obs scatter + RMSE/R²
python tools/plot_seasonal.py --stations 20 60 100 140
# → OUT/figures/fig_seasonal.png  +  OUT/tables/validation_seasonal.csv

# ---- 4. Hydrodynamics & transport diagnostics (3×5 panel grid) ----
python tools/plot_hydrodynamics.py                # model-only
python tools/plot_hydrodynamics.py --with-obs     # + obs + metrics
# → OUT/figures/fig_hydrodynamics.png  +  OUT/tables/validation_hydrodynamics.csv

# ---- 5. Reaction rates and limitation factors ----
python tools/plot_reactions.py                    # all rate groups
python tools/plot_reactions.py --mode limitations # nutrient limitation only
# → OUT/figures/fig_reactions.png
```

## Scripts

### Core

| Script | Purpose |
|--------|---------|
| `cgem_read.py` | I/O module — reads model output (CSV / CGEMBIN), variable metadata, season classification |
| `check_inputs.py` | Pre-run input data quality check + diagnostic figure |

### Plotting & Validation (4 scripts — one figure each)

| Script | Figure | Description | Obs support |
|--------|--------|-------------|:-----------:|
| `plot_profiles.py` | `fig_profiles.png` | Longitudinal WQ profiles: mean ± std, dry/wet overlays, obs scatter, RMSE/R² | `--with-obs` |
| `plot_seasonal.py` | `fig_seasonal.png` | Time-series grid (params × stations): Godin-filtered mean, min/max envelope, obs scatter | `--with-obs` |
| `plot_hydrodynamics.py` | `fig_hydrodynamics.png` | 3×5 grid: geometry, tidal range, salinity, particle trajectories, transport | `--with-obs` |
| `plot_reactions.py` | `fig_reactions.png` | Reaction rates by process group: mean ± std with dry/wet seasonal overlays | — |

### Utility

| Script | Purpose |
|--------|---------|
| `convert_obs_to_model_units.py` | Convert field observations to model-unit copies without overwriting the originals |

## Model Output Files

The model writes to `OUT/` with one file per variable. By default these are
CSV files; if `output_storage_format = BINARY` is set in `INPUT/params.txt`,
the same variables are written as `.bin` instead.

| Category | Files | Unit |
|----------|-------|------|
| **Water quality** | `Phy1.csv/.bin`, `Phy2.csv/.bin`, `Si.csv/.bin`, `NO3.csv/.bin`, `NH4.csv/.bin`, `PO4.csv/.bin`, `PIP.csv/.bin`, `O2.csv/.bin`, `TOC.csv/.bin`, `Sal.csv/.bin`, `SPM.csv/.bin` | mmol m⁻³ (or PSU, kg/m³) |
| **Carbonate** | `DIC.csv/.bin`, `AT.csv/.bin`, `pCO2.csv/.bin`, `PH.csv/.bin`, `CO2.csv/.bin` | mmol m⁻³, µeq/L, µatm |
| **Hydrodynamics** | `velocity.csv/.bin`, `waterDepth.csv/.bin`, `Discharge.csv/.bin`, `disp.csv/.bin`, `tau_b.csv/.bin` | m/s, m, m³/s |
| **Reaction rates** | `Reaction_NPP_NO3.csv/.bin`, `Reaction_nitrification.csv/.bin`, etc. | mmol m⁻³ s⁻¹ |
| **Limitation** | `Diag_fN.csv/.bin`, `Diag_fP.csv/.bin`, `Diag_fSi.csv/.bin`, `Diag_fI.csv/.bin`, `Diag_KD.csv/.bin` | dimensionless (0–1) |
| **Geometry** | `geometry.csv` | m, m² |

### CSV / Binary Format

```
Date,Cell_1,Cell_2,...,Cell_M
2017-01-01 00:00,val_1,val_2,...,val_M
```

Binary files use the CGEM per-variable stream format that is read directly by
`tools/cgem_read.py`; you do not need to convert them before plotting.

## Common Options

| Flag | Description |
|------|-------------|
| `--output-dir PATH` | Model output directory (default: `OUT`) |
| `--with-obs` | Enable obs overlay using default file in `INPUT/Validation/` |
| `--obs-file PATH` | Observation CSV file (overrides `--with-obs`) |
| `--vars VAR1 VAR2` | Variables to plot/validate |
| `--save FILE` | Save figure to file instead of displaying |
| `--save-dir PATH` | Directory for output figures/tables |

## Programmatic Use

```python
from tools.cgem_read import load_output, load_geometry, list_variables
from tools.cgem_read import classify_season, get_variable_info

# List available variables
print(list_variables())

# Load salinity output as DataFrame
sal = load_output("Sal")
print(sal.head())

# Classify dates into dry/wet seasons
seasons = classify_season(sal.index)

# Load geometry
geom = load_geometry()

# Get variable metadata
name, unit = get_variable_info("O2")  # ("Dissolved oxygen", "mmol O₂ m⁻³")
```

## Customizing Season Definitions

The default tropical season split is:
- **Dry:** November–April (months 11, 12, 1, 2, 3, 4)
- **Wet:** May–October (months 5, 6, 7, 8, 9, 10)

To change this for validation overlays, set `CGEM_DRY_MONTHS` before running
the tools (for example `11,12,1,2,3,4`). If an observation table already
contains a `Season` column, that explicit label is used instead.
