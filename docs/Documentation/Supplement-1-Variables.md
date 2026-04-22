# Academic Research Model Variables and Parameters

## Core State Variables for Estuarine Research Applications

| Variable           | Description                                         | Units     | Array Type   |
|--------------------|-----------------------------------------------------|-----------|--------------|
| `level[i]`         | Water surface elevation above datum                 | m         | Odd indices  |
| `velocity[i]`      | Depth-averaged velocity                             | m/s       | Even indices |
| `elevation[i]`     | Channel bed elevation relative to datum             | m         | All indices  |
| `waterDepth[i]`    | Total water depth (totalArea/width)                 | m         | All indices  |
| `totalArea[i]`     | Total cross-sectional area                          | m²        | All indices  |
| `freeArea[i]`      | Free-surface cross-sectional area                   | m²        | All indices  |
| `totalAreaOld[i]`  | Previous-timestep cross-sectional area              | m²        | All indices  |
| `width[i]`         | Channel width                                       | m         | All indices  |
| `riverbed_depth[i]`| Riverbed depth below datum                          | m         | All indices  |
| `rs[i]`            | Storage width ratio                                 | -         | All indices  |
| `Chezy[i]`         | Chezy roughness coefficient                         | m^(1/2)/s | All indices  |
| `FRIC[i]`          | Friction factor (1/Chezy[i]²)                       | -         | All indices  |

## Hydrodynamic Working Variables

| Variable                   | Description                               | Units | Notes                                                                              |
|----------------------------|-------------------------------------------|-------|------------------------------------------------------------------------------------|
| `C[i][j]`                  | Coefficient matrix for tridiagonal system | -     | j=1: lower diagonal, j=2: main diagonal, j=3: upper diagonal, j=4: right-hand side |
| `Z[i]`                     | Right-hand side of equations              | -     | Storage for tridiagonal system                                                     |
| `E[i]`, `Y[i]`             | Convergence test arrays                   | -     | For water levels and velocities                                                    |
| `slope[i]`                 | Channel bed slope                         | -     | Read from riverbed profile file                                                    |
| `tempFreeArea[i]`          | Temporary free-surface area               | m²    | Used during iterations                                                             |
| `tempVelocity[i]`          | Temporary velocity                        | m/s   | Used during iterations                                                             |
| `tidal_range_mean[i]`      | Mean tidal range                          | m     | Tracked per cell for diagnostics                                                   |
| `daily_min_level[i]`       | Daily minimum water level                 | m     | Tracked for tidal range computation                                                |
| `daily_max_level[i]`       | Daily maximum water level                 | m     | Tracked for tidal range computation                                                |
| `velocity_min[i]`          | Minimum velocity                          | m/s   | Tracked for flow reversal diagnostics                                              |
| `velocity_max[i]`          | Maximum velocity                          | m/s   | Tracked for flow reversal diagnostics                                              |

## Transport Variables

| Variable           | Description                                | Units   | Notes                                |
|--------------------|--------------------------------------------|---------|--------------------------------------|
| `v[s].c[i]`        | Concentration for variable s at location i | mass/m³ | Main state variable for each species |
| `v[s].clb`         | Lower (downstream) boundary concentration  | mass/m³ | Boundary condition                   |
| `v[s].cub`         | Upper (upstream) boundary concentration    | mass/m³ | Boundary condition                   |
| `v[s].avg[i]`      | Time-averaged concentration                | mass/m³ | For output and analysis              |
| `v[s].concflux[i]` | Cumulative mass flux                       | mass    | For mass balance accounting          |
| `v[s].advflux[i]`  | Advective flux component                   | mass/s  | For detailed flux analysis           |
| `v[s].disflux[i]`  | Dispersive flux component                  | mass/s  | For detailed flux analysis           |
| `disp[i]`          | Dispersion coefficient                     | m²/s    | Calculated in `Dispcoef()`           |
| `fl[i]`            | Advective flux at cell faces               | mass/s  | Calculated in `TVD()`                |
| `watflux[i]`       | Water flux                                 | m³/s    | For water balance                    |

## Biogeochemical Variables

| Variable       | Description                      | Units      | Chemical Representation    |
|----------------|----------------------------------|------------|----------------------------|
| `v[Phy1].c[i]` | Diatom biomass                   | mmol C/m³  | Phytoplankton group 1      |
| `v[Phy2].c[i]` | Non-diatom phytoplankton biomass | mmol C/m³  | Phytoplankton group 2      |
| `v[Si].c[i]`   | Dissolved silica                 | mmol Si/m³ | Si(OH)₄                    |
| `v[NO3].c[i]`  | Nitrate                          | mmol N/m³  | NO₃⁻                       |
| `v[NH4].c[i]`  | Ammonium                         | mmol N/m³  | NH₄⁺                       |
| `v[PO4].c[i]`  | Phosphate                        | mmol P/m³  | PO₄³⁻                      |
| `v[PIP].c[i]`  | Particulate inorganic phosphorus | mmol P/m³  | Adsorbed to SPM (auto-calc) |
| `v[O2].c[i]`   | Dissolved oxygen                 | mmol O₂/m³ | O₂                         |
| `v[TOC].c[i]`  | Total organic carbon             | mmol C/m³  | Various organic compounds  |
| `v[Sal].c[i]`  | Salinity                         | psu        | Dissolved salts            |
| `v[SPM].c[i]`  | Suspended particulate matter     | kg/m³      | Sediment (often plotted as g/m³ after scaling) |
| `v[DIC].c[i]`  | Dissolved inorganic carbon       | mmol C/m³  | CO₂ + HCO₃⁻ + CO₃²⁻        |
| `v[AT].c[i]`   | Total alkalinity                 | meq/m³     | Acid neutralizing capacity |
| `v[PH].c[i]`   | pH                               | -          | -log₁₀[H⁺]                 |
| `v[CO2].c[i]`


## Hydrodynamic Variables

| Parameter       | Unit | Suitable Range | Notes                                                  |
|-----------------|------|----------------|--------------------------------------------------------|
| **Water Level** | m    | -5 to 10       | Relative to datum                                      |
| **Velocity**    | m/s  | -5 to 5        | Negative indicates upstream flow                       |
| **Water Depth** | m    | 0.1 to 30      | Avoid depths < 0.1m to prevent numerical instabilities |
| **Width**       | m    | 10 to 2000     | Channel width                                          |
| **Discharge**   | m³/s | 1 to 50,000    | For river/estuary systems                              |

## State Variables (Concentrations)

| Parameter                      | Unit       | Suitable Range | Typical Freshwater | Typical Seawater | Notes                                                    |
|--------------------------------|------------|----------------|--------------------|------------------|----------------------------------------------------------|
| **Phytoplankton (Phy1, Phy2)** | mmol C/m³  | 0 to 500       | 5-50               | 1-30             | Diatoms and non-diatoms; >500 indicates bloom conditions |
| **Silica (Si)**                | mmol Si/m³ | 0 to 400       | 100-300            | 0-10             | Values >400 rare except in highly enriched systems       |
| **Nitrate (NO₃)**              | mmol N/m³  | 0 to 500       | 0-100              | 0-30             | Values >200 indicate high pollution                      |
| **Ammonium (NH₄)**             | mmol N/m³  | 0 to 200       | 0-10               | 0-1              | Values >50 indicate sewage/animal waste inputs           |
| **Phosphate (PO₄)**            | mmol P/m³  | 0 to 20        | 0-5                | 0-2              | Values >10 indicate high pollution                       |
| **Oxygen (O₂)**                | mmol O₂/m³ | 0 to 500       | 150-350            | 150-250          | 0-125 indicates hypoxic conditions; 250 mmol/m³ ≈ 8 mg/L |
| **Total Organic Carbon (TOC)** | mmol C/m³  | 0 to 2000      | 100-500            | 50-200           | Values >1000 indicate high organic loading               |
| **Salinity**                   | psu        | 0 to 35        | 0-0.5              | 30-35            | Critical for estuarine dynamics                          |
| **Suspended Matter (SPM)**     | kg/m³ (≈ g/L) | 0 to 1.0   | 0.005-0.1          | 0.001-0.01       | Public validation plots often display this as g/m³ (= mg/L) |

## Carbonate System Parameters

| Parameter                            | Unit      | Suitable Range | Typical Freshwater | Typical Seawater | Notes                                                               |
|--------------------------------------|-----------|----------------|--------------------|------------------|---------------------------------------------------------------------|
| **Dissolved Inorganic Carbon (DIC)** | mmol C/m³ | 100 to 2500    | 500-2000           | 1800-2200        | Depends on pH and alkalinity                                        |
| **Total Alkalinity (AT)**            | meq/m³    | 100 to 2500    | 500-2000           | 2200-2400        | Should be similar to DIC in most natural waters                     |
| **DIC/AT ratio**                     | -         | 0.9 to 1.1     | 0.9-1.05           | 0.95-1.05        | **Critical parameter**: values >1.1 or <0.9 often indicate problems |
| **pH**                               | -         | 6.5 to 9.0     | 6.8-8.2            | 7.8-8.3          | Log scale; crucial for aquatic life                                 |
| **CO₂**                              | mmol C/m³ | 5 to 100       | 10-50              | 5-15             | Dissolved CO₂ concentration                                         |
| **pCO₂**                             | μatm      | 100 to 5000    | 400-2000           | 200-600          | Partial pressure of CO₂; 400 μatm is atmospheric                    |

## Calculated Fluxes and Process Rates

| Parameter                    | Unit         | Suitable Range | Notes                                                         |
|------------------------------|--------------|----------------|---------------------------------------------------------------|
| **Primary Production (NPP)** | mmol C/m³/d  | 0 to 200       | Net primary production; >100 indicates high productivity      |
| **Respiration/Degradation**  | mmol C/m³/d  | 0 to 100       | Organic matter degradation; should balance with NPP over time |
| **Nitrification**            | mmol N/m³/d  | 0 to 20        | NH₄ to NO₃ conversion; reduces alkalinity                     |
| **Denitrification**          | mmol N/m³/d  | 0 to 10        | NO₃ reduction; increases alkalinity                           |
| **O₂ Air-Water Exchange**    | mmol O₂/m³/d | -100 to 100    | Negative values indicate O₂ uptake from atmosphere            |
| **CO₂ Air-Water Exchange**   | mmol C/m³/d  | -50 to 50      | Negative values indicate CO₂ uptake from atmosphere           |
| **Dispersion Coefficient**   | m²/s         | 1 to 1000      | Higher in estuaries with strong tidal influence               |

## Reaction-Specific Rates

| Process                     | Unit | Typical Range | Notes                              |
|-----------------------------|------|---------------|------------------------------------|
| **Phytoplankton Growth**    | /d   | 0.1 to 2.0    | Temperature dependent              |
| **Phytoplankton Mortality** | /d   | 0.01 to 0.2   | Temperature dependent              |
| **TOC Degradation**         | /d   | 0.001 to 0.1  | Depends on organic matter lability |
| **Nitrification Rate**      | /d   | 0.01 to 0.5   | Depends on temperature and O₂      |
| **Denitrification Rate**    | /d   | 0.001 to 0.1  | Depends on NO₃ and O₂ inhibition   |

## Optional spatially selective Phy1 loss parameters

These parameters are **optional** and default to **0 (disabled)** unless set in `INPUT/params.txt`.
They apply only to Phy1 (diatoms) as additional mortality terms intended to improve spatial localization
in reaches where growth controls alone can be insufficient.

| Parameter | Units | Meaning | Notes |
|---|---:|---|---|
| `kmort_sal_Phy1` | s⁻¹ / PSU^p | amplitude for salinity-stress mortality | Applied only when `Sal > sal_ref_Phy1` |
| `sal_ref_Phy1` | PSU | salinity threshold for stress | Set to brackish onset (site-specific) |
| `sal_power_Phy1` | – | exponent for salinity stress | Use 1–2 unless strong evidence |
| `kmort_tau_Phy1` | s⁻¹ / Pa^p | amplitude for shear-stress mortality | Optional; should be enabled only if high-shear coincides with bloom suppression |
| `tau_ref_Phy1` | Pa | bed-shear threshold | Site-specific |
| `tau_power_Phy1` | – | exponent for shear stress | Use 1–2 |
| `kmort_photo_Phy1` | s⁻¹ / (W m⁻²)^p | amplitude for depth-mean-light photostress | Optional; may be less selective in 1D models |
| `Imean_ref_Phy1` | W/m² | depth-mean irradiance threshold | Used only if `kmort_photo_Phy1>0` |
| `Imean_power_Phy1` | – | exponent for photostress | Use 1–2 |
| `kmort_clear_Phy1` | s⁻¹ / (m⁻¹)^p | amplitude for clarity (KD-based) stress | Often the most selective option in vertically integrated 1-D applications |
| `KD_ref_Phy1` | m⁻¹ | KD threshold below which water is “too clear” | Clear water ⇒ deeper light/UV penetration |
| `KD_power_Phy1` | – | exponent for clarity stress | Use 1–2 |
| `kmort_extra_max_Phy1` | s⁻¹ | cap on added mortality | 0 = no cap; useful for numerical safety |

## Critical Ratio Values for Validating Results

| Ratio       | Suitable Range | Warning Signs                                                          |
|-------------|----------------|------------------------------------------------------------------------|
| **DIC:AT**  | 0.9 to 1.1     | Values >1.1 indicate DIC accumulation or AT depletion                  |
| **O₂:DIC**  | 0.1 to 0.3     | Values far outside this range may indicate imbalance in metabolism     |
| **N:P**     | 10 to 30       | Based on Redfield ratio (~16); values outside range suggest limitation |
| **DIC:TOC** | 3 to 20        | Low values may indicate organic matter accumulation                    |
