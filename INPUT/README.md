# INPUT — Saigon River Example Dataset

This folder contains the complete forcing, boundary-condition, calibration, and
validation data for a two-year (2017–2018) simulation of the Saigon River
estuary (202 km, tropical monsoon climate).  It serves as a **working example**
that you can adapt to any tidal estuary.

## Folder Structure

```
INPUT/
├── params.txt                  Model parameters (calibrated)
├── config_input.txt            Boundary/tributary file mapping
├── Boundary/
│   ├── UB/                     Upstream boundary (discharge, salinity, WQ)
│   ├── UB_stable/              Stable upstream endmember BCs
│   ├── LB/                     Downstream boundary (tidal elevation, light, WQ)
│   └── wind.csv                Wind speed (daily)
├── Geometry/
│   └── river_bed_profile.txt   Longitudinal bed elevation
├── Tributaries/
│   ├── Dongnai/                Major tributary (river; cellIndex 31)
│   ├── Canals/                 Urban canal system (cellIndex 37)
│   ├── VamThuat/               Urban tributary (cellIndex 45)
│   └── ThiTinh/                Peri-urban tributary (cellIndex 61)
├── carbonate/                  Derived DIC/AT from field pCO₂ + alkalinity
├── Calibration/
│   ├── calibration_parameters.csv    Parameter bounds per calibration stage
│   ├── targets_seasonal_unified.csv  Seasonal WQ observations for scoring
│   └── targets_tidal.csv             Tidal range observations for Stage 1
└── Validation/
    └── *.csv / *.xlsx          Independent observation datasets
```

All input files — raw observations, derived boundary conditions, and carbonate
chemistry computed from field data — are contained in this single `INPUT/`
folder and referenced via `config_input.txt`.

## How to Adapt for Your Own Estuary

1. **Geometry** — Replace `river_bed_profile.txt` with your estuary's
  longitudinal bathymetry.  Update `EL` (length) and `DELXI` (grid spacing)
  in `params.txt`.

2. **Boundary conditions** — Prepare two-column numeric forcing files for
  upstream (discharge, salinity, water-quality concentrations) and downstream
  (tidal elevation, salinity, marine WQ).  Update file paths in  `config_input.txt`.

3. **Tributaries** — Add or remove tributary blocks in `config_input.txt`.
  Each tributary needs a `cellIndex` (grid cell where it enters) and CSV files
  for discharge and concentrations.  Set `numTributaries` accordingly.

4. **Parameters** — Adjust segment geometry (`B1`, `LC1`, `B2`, `LC2`, etc.)
  and tidal amplitude (`AMPL`).  Then re-calibrate using the staged approach:
  Stage 1 (hydrodynamics) → Stage 2 (sediment) → Stage 3 (biogeochemistry) →
  Stage 4 (carbonate).

5. **Validation** — Place independent field observations in `Validation/`.

## Forcing file format

The model input parser reads **two-column numeric time series**:

```
1,123.4
2,125.1
3,124.8
```

- The first column is a numeric time index, not a calendar-date string.
- The bundled example files use sequential counters (`1, 2, 3, ...`).
- Resolution is defined by the relevant block in `config_input.txt`
  (`resolution=daily` or `resolution=hourly`).
- Daily files must cover the full simulation period; hourly files (for example
  elevation and light) must cover the full sub-daily record used by the run.

If your source observations are timestamped by date, preprocess them into this
numeric two-column format before pointing `config_input.txt` to the file.

## Units Convention

| Variable    | Unit       | Notes                       |
| ----------- | ---------- | --------------------------- |
| Discharge   | m³/s       | Daily mean                  |
| Salinity    | PSU        |                             |
| Phy1/Phy2   | mmol C/m³  | Diatom / non-diatom biomass |
| NO₃, NH₄    | mmol N/m³  |                             |
| PO₄         | mmol P/m³  |                             |
| Si          | mmol Si/m³ |                             |
| O₂          | mmol O₂/m³ |                             |
| TOC         | mmol C/m³  |                             |
| SPM         | g/m³       |                             |
| DIC         | mmol C/m³  |                             |
| AT          | mmol/m³    | Total alkalinity            |
| Temperature | °C         |                             |
| Wind        | m/s        |                             |
| Light       | W/m²       | PAR at surface              |
| Elevation   | m          | Tidal water level           |

