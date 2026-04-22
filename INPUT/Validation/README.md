# Validation Data — C-GEM

This directory holds observation data for validating C-GEM model outputs.
The validation tools in `tools/` work **without any observation files**
(model-only visualisation).  When observation files are provided via the
`--obs-file` command-line argument, the tools overlay observed data and
compute validation metrics.

The repository already ships three ready-to-use observation tables:

- `obs_longitudinal_profiles.csv`
- `obs_timeseries.csv`
- `obs_tidal_range.csv`

These bundled files are already in **model units** and can be used directly
with the `--with-obs` flag.

> **Key rule:** All observation values must be provided in **model units**
> (see table below).  The tools perform **no unit conversion**.

---

## Model Units Reference

| Variable | Model name | Model unit | Notes |
|----------|-----------|------------|-------|
| Salinity | `Sal` | PSU | — |
| Dissolved oxygen | `O2` | mmol O₂ m⁻³ | |
| Ammonium | `NH4` | mmol N m⁻³ | |
| Nitrate | `NO3` | mmol N m⁻³ | |
| Phosphate | `PO4` | mmol P m⁻³ | |
| Total organic carbon | `TOC` | mmol C m⁻³ | |
| Dissolved inorganic carbon | `DIC` | mmol C m⁻³ | |
| Dissolved CO₂ | `CO2` | mmol C m⁻³ | |
| Diatoms | `Phy1` | mmol C m⁻³ | Phytoplankton biomass as carbon |
| Non-siliceous phyto | `Phy2` | mmol C m⁻³ | |
| Dissolved silica | `Si` | mmol Si m⁻³ | |
| **Suspended matter** | `SPM` | **kg m⁻³** | Not mg/L — multiply field mg/L × 0.001 |
| **Total alkalinity** | `AT` | **µeq L⁻¹** | Micro-equivalents per litre |
| pCO₂ | `pCO2` | µatm | |
| pH | `PH` | — | Dimensionless |

### Common field-to-model conversions

If your observations are in field units (mg/L, µg/L, etc.), convert them
before placing them in the CSV:

| From (field) | To (model) | Formula |
|-------------|-----------|---------|
| mg O₂/L → mmol O₂/m³ | O2 | value / 0.032 |
| mg N/L → mmol N/m³ | NH4, NO3 | value / 0.014007 |
| mg P/L → mmol P/m³ | PO4 | value / 0.030974 |
| mg C/L → mmol C/m³ | TOC, DIC | value / 0.012011 |
| mg Si/L → mmol Si/m³ | Si | value / 0.028085 |
| mg/L → kg/m³ | SPM | value × 0.001 |
| µg Chl-a/L → mmol C/m³ | Phy1 | value × C:Chl ratio / MW_C |
| mg CaCO₃/L → µeq/L | AT | value / 0.050045 |

---

## File Formats

### 1. Longitudinal profiles (`--obs-file` for `plot_profiles.py`)

```csv
Distance_km, Season, Sal, O2, NH4, NO3, PO4, TOC, SPM, Phy1, Si, DIC, AT, pCO2, PH
0,           dry,    25,  280, 5,  10,  1.5, 200, 0.05, 2,   50, 1800, 2400, 500, 7.8
10,          wet,    18,  250, 8,  15,  2.0, 300, 0.08, 3,   80, 2000, 2200, 800, 7.5
```

- **Column names** must match the model variable names exactly.
- Only include columns you have data for — missing columns are ignored.
- `Distance_km`: distance from estuary mouth (0 = mouth, increasing upstream).
- `Season`: `dry` or `wet` (optional; if absent, all data is treated together).

### 2. Time series (`--obs-file` for `plot_seasonal.py`)

```csv
Date, Station, Distance_km, Sal, O2, NH4, NO3
2017-03-15, Station_A, 40, 22.5, 280, 5.2, 12
2017-04-20, Station_A, 40, 20.1, 300, 4.8, 10
```

- `Date`: sampling date (YYYY-MM-DD).
- `Station`: station identifier (used for grouping in figures).
- `Distance_km`: station location.
- Variable columns: values in model units.

### 3. Tidal range (`--obs-file` for `plot_hydrodynamics.py`)

```csv
Distance_km, Tidal_range_m, Salinity, Season
0,           3.2,           25.0,     dry
10,          2.8,           20.1,     wet
```

- `Tidal_range_m`: observed tidal range in metres.
- `Salinity`: salinity in PSU (optional).
- `Season`: `dry` or `wet` (optional).

---

## Distance Convention

`Distance_km` is measured from the **estuary mouth** (0 km) increasing
**upstream**.  This matches the model grid where Cell 1 = 0 km (mouth)
and Cell N = upstream freshwater boundary.

---

## Season Convention

| Season | Months |
|--------|--------|
| **dry** | November–April (11, 12, 1, 2, 3, 4) |
| **wet** | May–October (5, 6, 7, 8, 9, 10) |

Adjust to your local conditions by modifying `DRY_MONTHS` / `WET_MONTHS`
with the `CGEM_DRY_MONTHS` environment variable, or by providing an explicit
`Season` column in the observation table.

---

## Preparing Validation Data for a New Estuary

1. Collect observations at monitoring stations along the estuary.
2. Convert all values to **model units** (see table above).
3. Assign `Distance_km` to each station (km from mouth, increasing upstream).
4. Assign `Season` using your local dry/wet convention.
5. Use exact **model variable names** as column headers.
6. Leave cells empty for missing values (do not write "NA" or -999).
7. Save as CSV (UTF-8 encoding, comma separator) in this folder.

---

## Running the Tools

```bash
# Model-only (no observations needed):
python tools/plot_profiles.py --output-dir OUT
python tools/plot_seasonal.py --output-dir OUT --vars Sal O2
python tools/plot_hydrodynamics.py --output-dir OUT

# With the bundled observation tables in INPUT/Validation/:
python tools/plot_profiles.py --with-obs
python tools/plot_seasonal.py --with-obs --vars Sal O2
python tools/plot_hydrodynamics.py --with-obs

# With a custom observation file:
python tools/plot_profiles.py --obs-file path/to/my_profiles.csv
python tools/plot_seasonal.py --obs-file path/to/my_timeseries.csv
python tools/plot_hydrodynamics.py --obs-file path/to/my_hydro.csv

# Convert field-unit tables to model-unit copies (optional helper):
python tools/convert_obs_to_model_units.py
```

Output figures are saved to `OUT/figures/` and metric tables to `OUT/tables/`.
