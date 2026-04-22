# Biogeochemical Module

The biogeochemical module implements a 17-species reaction network coupling
nitrogen, carbon, phosphorus, silica, and oxygen cycles. The framework builds
on Volta et al. (2014, 2016) with tropical adaptations for
temperature-dependent kinetics, sediment-phosphorus interactions, and CO₂ gas
exchange.

## Species Network

```c
typedef enum {
    Phy1,   // Diatoms [mmol C/m³]
    Phy2,   // Non-diatoms [mmol C/m³]
    Si,     // Dissolved silica [mmol Si/m³]
    NO3,    // Nitrate [mmol N/m³]
    NH4,    // Ammonium [mmol N/m³]
    PO4,    // Phosphate [mmol P/m³]
    PIP,    // Particulate inorg. phosphorus [mmol P/m³] — diagnostic
    O2,     // Dissolved oxygen [mmol O₂/m³]
    TOC,    // Total organic carbon [mmol C/m³]
    Sal,    // Salinity [PSU]
    SPM,    // Suspended particulate matter [kg/m³]
    DIC,    // Dissolved inorganic carbon [mmol C/m³]
    AT,     // Total alkalinity [µeq/L]
    pCO2,   // Partial pressure CO₂ [µatm] — diagnostic
    PH,     // pH [-] — diagnostic
    CO2,    // Dissolved CO₂ [mmol C/m³] — diagnostic
    CHEM_COUNT  // = 17
} Chem;
```

Four species (PIP, pCO₂, pH, CO₂) are computed diagnostically at each timestep
— they are not independently transported.

## Reaction-Transport Equation

Each transported species $C_i$ evolves as:

$$\frac{\partial C_i}{\partial t} + U\frac{\partial C_i}{\partial x} = \frac{\partial}{\partial x}\!\left(D\frac{\partial C_i}{\partial x}\right) + R_i + S_i$$

where $R_i$ is the biogeochemical source/sink and $S_i$ represents external inputs. All rates are temperature-dependent:

$$k(T) = k(T_\text{ref}) \times Q_{10}^{(T - T_\text{ref})/10}$$

## Primary Production

Two phytoplankton groups (Phy1 = diatoms, Phy2 = non-diatoms) grow via Monod kinetics with multiplicative limitation:

$$\mu = \mu_\text{max} \times f(T) \times f(I) \times \min\!\bigl(f_N,\; f_P,\; [f_{Si}]\bigr)$$

### Light limitation

$$f(I) = \frac{\bar{I}}{K_I + \bar{I}}, \qquad \bar{I} = \frac{I_0}{K_d \cdot h}\bigl(1 - e^{-K_d \cdot h}\bigr)$$

Light attenuation includes background, SPM, CDOM, and self-shading:

$$K_d = k_\text{bg} + k_\text{spm} \cdot \text{SPM} + k_\text{CDOM} \cdot \text{TOC} + k_\text{Phy} \cdot (\text{Phy1} + \text{Phy2})$$

### Nutrient limitation (Liebig's law of the minimum)

| Factor | Expression |
|--------|-----------|
| $f_N$ | $\dfrac{[\text{NH}_4] + [\text{NO}_3]}{K_N + [\text{NH}_4] + [\text{NO}_3]}$ |
| $f_P$ | $\dfrac{[\text{PO}_4]}{K_P + [\text{PO}_4]}$ |
| $f_{Si}$ (diatoms only) | $\dfrac{[\text{Si}]}{K_{Si} + [\text{Si}]}$ |

### NH₄ vs NO₃ preference

Phytoplankton preferentially uptake NH₄. NO₃ uptake is inhibited when NH₄ is abundant:

$$f_\text{NH4-pref} = \frac{[\text{NH}_4]}{K_\text{NH4} + [\text{NH}_4]}, \qquad f_\text{NO3} = \frac{[\text{NO}_3]}{K_\text{NO3} + [\text{NO}_3]} \cdot \frac{K_\text{NH4}}{K_\text{NH4} + [\text{NH}_4]}$$

### Tropical temperature cap

A Gaussian penalty is applied above 32 °C to represent heat stress (`tropical_temperature_limit()` in `biogeout.c`).

## Phytoplankton Losses

| Process | Rate expression |
|---------|----------------|
| Linear mortality | $k_\text{mort} \cdot \text{Phy}$ |
| Quadratic mortality | $k_\text{mort2} \cdot \text{Phy}^2$ |
| Benthic grazing | $k_\text{benth} \cdot \text{Phy}$ (spatially selective) |

Dead phytoplankton biomass enters the TOC pool.

## Nitrogen Cycle

| Process | Equation | O₂ dependency |
|---------|----------|---------------|
| Nitrification | $R_\text{nit} = k_\text{nit} \cdot [\text{NH}_4] \cdot \dfrac{[\text{O}_2]}{K_\text{O2} + [\text{O}_2]}$ | Aerobic |
| Denitrification | $R_\text{den} = k_\text{den} \cdot [\text{NO}_3] \cdot \dfrac{K_\text{inO2}}{K_\text{inO2} + [\text{O}_2]}$ | Anaerobic |
| Ammonification | $R_\text{amm} = k_\text{ox} \cdot [\text{TOC}] \cdot \dfrac{[\text{O}_2]}{K_\text{O2} + [\text{O}_2}} \cdot r_{N:C}$ | Aerobic |

Stoichiometric coupling uses Redfield ratios (C:N:P = 106:16:1 by atoms).

## Carbon Cycle

### Organic carbon

TOC is produced by phytoplankton mortality and consumed by aerobic degradation and CBOD fast oxidation:

$$R_\text{TOC} = \text{mortality sources} - k_\text{ox} \cdot [\text{TOC}] \cdot f(\text{O}_2) - k_\text{cbod} \cdot [\text{TOC}]$$

CDOM photo-degradation ($k_\text{CDOM}$) provides an additional first-order light-dependent DOC loss.

### Dissolved inorganic carbon

DIC changes from primary production (uptake), respiration (release), nitrification, and CO₂ gas exchange.

## Oxygen Dynamics

$$\frac{d[\text{O}_2]}{dt} = \text{NPP} \cdot r_{\text{O}_2:\text{C}} - \text{nitrification} \cdot r_{\text{O}_2:\text{N}} - \text{respiration} - \frac{\text{SOD}}{h} + F_{\text{O}_2}$$

where SOD is benthic sediment oxygen demand and $F_{\text{O}_2}$ is air-water gas exchange.

### Gas exchange

$$F = k_\text{gas} \cdot (C_\text{sat} - C_\text{water})$$

Gas transfer velocity uses a wind-speed parameterisation scaled by Schmidt number:

$$k_\text{gas} = k_{600} \cdot \left(\frac{Sc}{600}\right)^{-0.5}, \qquad k_{600} = a_1 + a_2 U_{10}^2$$

## Phosphorus Cycle

Dissolved PO₄ is taken up by phytoplankton and released by organic matter mineralisation. PIP (particulate inorganic phosphorus) is computed via Langmuir equilibrium with SPM:

$$\text{PIP}_\text{eq} = P_\text{ac} \cdot \text{SPM} \cdot \frac{[\text{PO}_4]}{[\text{PO}_4] + K_{ps}}$$

Relaxation toward equilibrium is implicit:

$$\text{PIP}^{n+1} = \frac{\text{PIP}^n + k_\text{ads}\,\Delta t \cdot \text{PIP}_\text{eq}}{1 + k_\text{ads}\,\Delta t}$$

## Carbonate Chemistry

DIC and AT are transported; pCO₂, pH, and CO₂ are computed diagnostically each timestep by solving the carbonate equilibrium system:

$$[\text{HCO}_3^-] = \frac{K_1 [\text{CO}_2]}{[\text{H}^+]}, \qquad [\text{CO}_3^{2-}] = \frac{K_1 K_2 [\text{CO}_2]}{[\text{H}^+]^2}$$

$$\text{AT} = [\text{HCO}_3^-] + 2[\text{CO}_3^{2-}] + [\text{OH}^-] - [\text{H}^+]$$

Given DIC and AT, pH is found iteratively. The dissociation constants $K_1$, $K_2$, $K_w$ are temperature- and salinity-dependent (Millero, 1995). CO₂ air-water exchange uses:

$$F_{\text{CO}_2} = k_\text{gas} \cdot K_H \cdot (p\text{CO}_{2,\text{water}} - p\text{CO}_{2,\text{atm}})$$

where $K_H$ is Henry's constant and $p\text{CO}_{2,\text{atm}}$ defaults to 415 µatm.

## Reaction Rate Outputs

When `enable_reaction_output = 1`, reaction rates are written to
`OUT/Reaction_*` and `OUT/Diag_*` using the storage format selected by
`output_storage_format` (`.csv` or `.bin`):

| Output file | Content |
|-------------|---------|
| `Reaction_NPP_NO3` | Net primary production via NO₃ uptake |
| `Reaction_NPP_NH4` | Net primary production via NH₄ uptake |
| `Reaction_phydeath` | Phytoplankton mortality rate |
| `Reaction_adegradation` | Aerobic organic matter degradation |
| `Reaction_nitrification` | NH₄ → NO₃ oxidation rate |
| `Reaction_denitrification` | NO₃ reduction rate |
| `Reaction_O2_exchange` | Air-water O₂ flux |
| `Reaction_CO2_exchange` | Air-water CO₂ flux |
| `Reaction_DIC` | Net DIC source/sink |
| `Reaction_TA` | Net alkalinity source/sink |
| `Diag_fN`, `Diag_fP`, `Diag_fSi`, `Diag_fI` | Nutrient and light limitation factors |
| `Diag_KD` | Light attenuation coefficient |

## Code Structure

| Function | Purpose |
|----------|---------|
| `Biogeo(t)` | Main entry: loops over cells |
| `computePrimaryProduction()` | Growth rates + limitation factors |
| `computeBiogeochemicalReactions()` | N, C, P, O₂ reaction rates |
| `computeCarbonateChemistry()` | pH solver, pCO₂, CO₂ |
| `updateBiogeochemicalState()` | Apply rates to concentrations |
| `Rates(co, s, t)` | Write reaction diagnostics |

All implemented in `biogeo.c` (~1200 lines).

## Limitations

- 1D vertically mixed: no stratification or vertical gradients
- Fixed Redfield stoichiometry (C:N:P = 106:16:1)
- No sulfur, iron, or trace metal cycles
- Single gas-transfer parameterisation (no wave/fetch effects)
- Equilibrium carbonate chemistry (no kinetic disequilibrium)
- Lumped microbial kinetics (no explicit microbial community dynamics)

## References

- Volta, C., et al. (2014). *Biogeochemistry of the Scheldt estuary.* Biogeosciences, 11, 5007-5025.
- Volta, C., et al. (2016). *CO₂ dynamics in the Scheldt estuary.* J. Marine Systems, 157, 62-76.
- Regnier, P., et al. (1997). *Modelling estuarine biogeochemistry.* Aquatic Geochemistry, 3, 77-116.
- Millero, F.J. (1995). *Thermodynamics of the CO₂ system in the oceans.* Geochimica et Cosmochimica Acta, 59, 661-677.
