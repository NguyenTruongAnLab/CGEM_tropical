# Introduction to C-GEM

C-GEM is a 1D reactive-transport model for simulating biogeochemical processes in tidal estuaries. It couples hydrodynamics (Saint-Venant equations), solute transport (TVD advection + dispersion), and a 16-species biogeochemical network covering nitrogen, carbon, phosphorus, oxygen, and carbonate chemistry.

## What C-GEM Does

- Solves tidal hydrodynamics on a staggered grid with a semi-implicit scheme
- Transports 16 chemical species using TVD advection + Crank-Nicolson dispersion
- Computes biogeochemical transformations: primary production, nitrification, denitrification, organic matter degradation, air-water gas exchange, and carbonate equilibrium
- Supports tributary inputs, sediment dynamics (erosion/deposition), and phosphorus adsorption
- Runs seasonal-to-multi-year simulations on a single CPU core in minutes

## Design Philosophy

C-GEM uses physics-based simplifications to reduce data requirements:

| Principle | Implementation |
|-----------|---------------|
| Exponential geometry | Width follows $B(x) = B_0 e^{-x/L_b}$ (Savenije, 2005) |
| Predictive dispersion | Van der Burgh coefficient relates dispersion to tidal mixing |
| Literature-based kinetics | Biogeochemical rates from Volta et al. (2014, 2016) compilation |
| Minimal input data | Operates with basic geometry, tidal forcing, and river discharge |

This makes C-GEM practical for estuaries where detailed bathymetry or continuous monitoring data are unavailable.

## When to Use C-GEM (and When Not To)

**Good fit:**

- Process-oriented studies of estuarine biogeochemistry
- Seasonal water-quality assessment
- Carbon flux and CO₂ exchange estimation
- Scenario analysis (nutrient loading, climate change)
- Teaching and academic research

**Consider alternatives when you need:**

| Need | Alternative |
|------|-------------|
| Lateral (2D/3D) resolution | Delft3D, TELEMAC-2D |
| Vertical stratification | MIKE 3, ROMS |
| Complex delta networks | Delft3D-FM, MIKE 11 |
| Detailed morphodynamics | Delft3D, TELEMAC |

## Development History

| Period | Contributors | Scope |
|--------|-------------|-------|
| 2013–2016 | Chiara Volta (ULB, Belgium), supervised by P. Regnier, S. Arndt, H.H.G. Savenije, G.G. Laruelle | Original framework for European estuaries |
| 2018–2021 | Nguyen Truong An (HCMUT, Vietnam), supervised by J. Nemery, N. Gratiot, J. Garnier, G.G. Laruelle | Tropical adaptation, urban discharge integration |
| 2025 | Modular rewrite | Automated calibration, improved numerics, portable C99 code |

## Architecture Overview

```
Main loop (main.c):
  for each timestep:
    Hyd(t)           → Saint-Venant hydrodynamics
    bgboundary(t)    → Boundary condition interpolation
    Transport(t)     → TVD advection + Crank-Nicolson dispersion
    Biogeo(t)        → Biogeochemical reactions
    updateSPM(t)     → Sediment erosion/deposition
    Output(t)        → Write results to OUT/
```

The modular source files map directly to these steps:

| Module | File | Purpose |
|--------|------|---------|
| Hydrodynamics | `hydrodynamics.c` | Saint-Venant solver |
| Transport | `transport.c` | TVD + Crank-Nicolson |
| Biogeochemistry | `biogeo.c` | Reaction network |
| Forcing | `bcforcing.c`, `biogeout.c` | Boundary conditions, temperature/light |
| I/O | `io_input.c`, `io_output.c` | File reading and output |
| Diagnostics | `diagnostics.c` | Runtime sanity checks |
| Calibration | `calibration.c`, `calibration_helpers.c` | NLopt-based optimization |

## Biogeochemical Species

C-GEM tracks 16 species across 6 categories:

| Category | Species | Role |
|----------|---------|------|
| Primary producers | Phy1 (diatoms), Phy2 (non-diatoms) | Phytoplankton biomass |
| Nutrients | NO₃, NH₄, PO₄, Si, PIP | Nutrient cycling |
| Oxygen | O₂ | Dissolved oxygen dynamics |
| Organic matter | TOC | Carbon cycling |
| Carbonate | DIC, AT, pCO₂, pH, CO₂ | Carbonate equilibrium + CO₂ exchange |
| Physical | Sal, SPM | Salinity and suspended sediment |

Four species (PIP, pCO₂, pH, CO₂) are computed diagnostically (not independently transported).

## Calibration Workflow

C-GEM uses a 4-stage automated calibration with NLopt:

| Stage | Targets | Parameters tuned |
|-------|---------|-----------------|
| 1 | Tidal range + salinity | Chézy, Rs, dispersion (D₀, K_VDB) |
| 2 | SPM profiles | Erosion/deposition rates |
| 3 | Nutrients, O₂, Chl-a | Biogeochemical rates |
| 4 | pCO₂, pH, DIC, AT | Gas exchange scaling |

See [Calibration Preparation Guide](Calibration_Preparation_Guide.md) for details.

## Quick Start

See [Getting Started](Getting-Started.md) for build instructions, running the model, and adapting to your own estuary.

## References

- Volta, C., et al. (2014). *Biogeochemistry of the Scheldt estuary: a model study.* Biogeosciences, 11, 5007-5025.
- Volta, C., et al. (2016). *CO₂ dynamics in the Scheldt estuary.* J. Marine Systems, 157, 62-76.
- Savenije, H.H.G. (2005). *Salinity and Tides in Alluvial Estuaries.* Elsevier.
- Gisen, J.I.A., et al. (2015). *Predictive dispersion in alluvial estuaries.* HESS, 19, 2495-2514.
