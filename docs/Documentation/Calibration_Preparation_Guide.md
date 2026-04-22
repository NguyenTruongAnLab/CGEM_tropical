# 📊 C-GEM Calibration Target Preparation Guide (Canonical CSVs)

This guide matches the **simplified calibration contract** implemented in `src/calibration.c`:

- Calibration is **stage-based**: `calibration_stage ∈ {1,2,3,4}`.
- Targets must be **canonical fixed-schema CSVs** and already in **model/scoring units**.
- No header permutations, no date parsing, and no unit conversion inside C calibration.

## Required inputs

### 1) Tidal targets (tidal range)

Configured by `calibration_tidal_targets_file` in `INPUT/params.txt`.

Accepted header (strict):

- `location_km,range_m`

Example:

- `location_km,range_m`
- `0.0,1.6`
- `25.0,1.2`
- `50.0,0.8`

### 2) Unified seasonal targets (all stages)

Configured by `calibration_wq_seasonal_targets_file`.

This is a **single long-form file**. Each row indicates the calibration stage via the `stage` column; during a run, only rows matching `calibration_stage` are scored.

Schema (strict):

`target_step_seconds,cell_i0,cell_i1,alpha,stage,variable,target_value,scale,weight`

- `target_step_seconds` is the model time coordinate (seconds). Targets are sampled at the first model step with $t \ge$ `target_step_seconds`.
- `cell_i0`, `cell_i1` are 1-based model cell indices.
- `alpha ∈ [0,1]` mixes cell endpoints: $model=(1-\alpha)\,v[i0] + \alpha\,v[i1]$.
- `stage ∈ {1,2,3,4}` selects which calibration stage this row belongs to.
- `variable` is a string; supported names include: `Sal`, `SPM`, `PO4`, `NH4`, `NO3`, `O2`, `TOC`, `pCO2`, and phytoplankton total via `Phy` / `Chla` / `Chl-a` / `ChlA`.
- `target_value` is the observed/scoring value in model units.
- `scale` is a per-variable scaling factor used in the loss (e.g., a typical magnitude).
- `weight` is the contribution weight (dimensionless).

### 3) Calibration parameter registry

Configured by `calibration_parameters_file`.

This file declares, per stage, the parameter names, bounds, and transforms used by the optimizer.

## Configuration knobs

In `INPUT/params.txt`:

- `calibration_tidal_targets_file = "..."`
- `calibration_stage = 1..4`
- `calibration_parameters_file = "..."`
- `calibration_wq_seasonal_targets_file = "..."`

## Notes

- Calibration is intentionally strict and configuration-driven for paper-ready reproducibility:
  - No legacy headers
  - No duplicate tidal locations
  - No environment-variable overrides for target parsing or unit conversion

- Scenario reproducibility (peer review): the model can optionally read a non-default parameter file
  by setting `CGEM_PARAMS_PATH` (defaults to `INPUT/params.txt`). This enables clean sensitivity runs
  without editing tracked input files.
- For Stage 1 calibration, both tidal targets and stage-1 rows in the unified seasonal file are typically used.

## Calibration stages — parameter details

### Stage 1: Hydrodynamics + Transport (9 parameters)

| # | Parameter | Description |
|---|-----------|-------------|
| 1 | `Chezy1` | Chézy roughness, segment 1 [m^½/s] |
| 2 | `Chezy2` | Chézy roughness, segment 2 [m^½/s] |
| 3 | `LC1` | Convergence length, segment 1 [m] |
| 4 | `LC2` | Convergence length, segment 2 [m] |
| 5 | `Rs1` | Storage width ratio, segment 1 [-] |
| 6 | `Rs2` | Storage width ratio, segment 2 [-] |
| 7 | `C_VDB` | Van der Burgh coefficient [-] |
| 8 | `D0_CORR` | Dispersion at mouth correction [m²/s] |
| 9 | `AMPL_CORR` | Tidal amplitude correction [m] |

Targets: tidal range (from tidal targets file) and salinity (stage-1 rows in unified file).

### Stage 2: Sediment + Phosphorus adsorption (15 parameters)

| # | Parameter | Description |
|---|-----------|-------------|
| 1–3 | `Mero1`, `tau_ero1`, `tau_dep1` | Erosion rate, critical shear (erosion/deposition), segment 1 |
| 4–6 | `Mero2`, `tau_ero2`, `tau_dep2` | Same, segment 2 |
| 7–9 | `Mero3`, `tau_ero3`, `tau_dep3` | Same, segment 3 |
| 10–12 | `Mero4`, `tau_ero4`, `tau_dep4` | Same, segment 4 |
| 13 | `P_ac` | Max P adsorption capacity [mg P / g SPM] |
| 14 | `K_ps` | Langmuir half-saturation for PO₄ [mg P/L] |
| 15 | `k_ads` | Adsorption rate constant [s⁻¹] |

Targets: SPM profiles (stage-2 rows).

### Stage 3: Eutrophication (20 parameters)

| # | Parameter | Description |
|---|-----------|-------------|
| 1 | `SCALE_KOX` | Aerobic degradation rate multiplier |
| 2 | `SCALE_KNIT` | Nitrification rate multiplier |
| 3 | `SCALE_KDENIT` | Denitrification rate multiplier |
| 4 | `SCALE_PB_PHY` | Max production rate multiplier (both Phy groups) |
| 5 | `SCALE_KMORT_PHY` | Mortality rate multiplier (both Phy groups) |
| 6 | `SCALE_KD_PHY` | Phy light attenuation multiplier |
| 7 | `SCALE_K_NH4_SWITCH` | NH₄ preference switch multiplier |
| 8 | `SCALE_KBG` | Background light attenuation multiplier |
| 9 | `SCALE_KSPM` | SPM light attenuation multiplier |
| 10 | `SCALE_WS_PHY` | Phytoplankton settling velocity multiplier |
| 11 | `SCALE_ALPHA` | Light efficiency (αPAR) multiplier |
| 12 | `SCALE_KN` | NO₃ half-saturation multiplier |
| 13 | `SCALE_KPO4` | PO₄ half-saturation multiplier |
| 14 | `SCALE_KSI` | Si half-saturation multiplier |
| 15 | `SCALE_KTOC` | TOC half-saturation multiplier |
| 16 | `SOD_RATE` | Sediment oxygen demand [mmol O₂/m²/s] |
| 17 | `KMORT2_PHY1` | Quadratic mortality, Phy1 [s⁻¹/(mmolC/m³)] |
| 18 | `KCBOD_FAST` | Fast labile-C degradation rate [s⁻¹] |
| 19 | `KCDOM` | CDOM light absorption [m⁻¹/(mgC/L)] |
| 20 | `CN_TOC` | Effective C:N for TOC mineralisation [molC/molN] |

Targets: Chl-a, O₂, NH₄, NO₃, PO₄, TOC profiles (stage-3 rows).

### Stage 4: Carbonate (1 parameter)

| # | Parameter | Description |
|---|-----------|-------------|
| 1 | `piston_velocity_scale` | Piston velocity multiplier [-] |

Targets: pCO₂, pH, DIC, AT (stage-4 rows).
