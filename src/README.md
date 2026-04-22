# C-GEM Tropical — Source Code

One-dimensional coupled hydrodynamic–biogeochemical model for data-sparse tropical estuaries, built on the C-GEM framework (Volta et al., 2014; Regnier et al., 2013). The model is configuration-driven: all site-specific information (geometry, tributaries, forcing data, parameters) lives in the `INPUT/` directory. The source code itself is generic.

## Compilation

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

The executable is placed in `bin/Release/C_GEM_Daily` (or `.exe` on Windows).

Optional build flags:
- `-DCGEM_ENABLE_CALIBRATION=ON` — link NLopt and enable automated calibration

## Source File Guide

### Core Physics

| File | Description | Key References |
|------|-------------|----------------|
| `hydrodynamics.c` | 1D Saint-Venant equations on a staggered grid; implicit friction, Chézy bottom shear, radiation upstream BC | Savenije (2012), Blayo & Debreu (2005), Orlanski (1976) |
| `transport.c` | Advection (TVD with CFL sub-stepping) + Crank–Nicolson dispersion; Savenije (2012) recursive dispersion with Van der Burgh coefficient; Fischer (1979) tributary mass injection | Savenije (2012), Fischer et al. (1979), Rutherford (1994), Gisen et al. (2015) |
| `bcforcing.c` | Boundary condition interpolation: tidal elevation (M2/S2/M4 or CSV), upstream discharge, wind, tributary discharge; warmup- forcing alignment | Savenije (2012) |
| `biogeo.c` | Biogeochemical reactions: NPZD (two phytoplankton groups), carbonate system, sediment dynamics, nutrient cycling, gas exchange | Jassby & Platt (1976), O'Connor & Dobbins (1958), Wanninkhof (1992), Fennel et al. (2006), Arndt et al. (2011), DOE (1994), Dickson (2007), Partheniades (1965), Krone (1962) |

### Initialisation and I/O

| File | Description |
|------|-------------|
| `main.c` | Simulation loop: warmup → analysis period; environment-variable overrides for run modes (`CGEM_TRANSPORT_ONLY`, `CGEM_WARMUP_DAYS`, etc.) |
| `init.c` / `init.h` | Parameter parsing (`params.txt`), geometry setup, boundary-condition file loading, initial-condition seeding |
| `io_input.c` | CSV/text file readers for forcing data, tributary concentrations, upstream chemistry |
| `io_output.c` | Daily and sub-daily output writers (per-variable CSV/BINARY, geometry snapshot, optional export hooks) |
| `biogeout.c` | Biogeochemical diagnostic outputs (reaction rates, limitation factors) |

### State and Configuration

| File | Description |
|------|-------------|
| `define.h` | Compile-time constants (`MAXM`, `G`), physical constants, enumerations |
| `variables.h` / `variables.c` | Global state variables, parameter declarations, default values |
| `file.h` | File-path helpers |

### Diagnostics and Utilities

| File | Description |
|------|-------------|
| `diagnostics.c` / `diagnostics.h` | Runtime validation: tracer non-negativity enforcement (clamp ≥ −1 × 10⁻⁴ to zero, abort below), numerical stability checks, mass-balance diagnostics |
| `utilities.c` / `utilities.h` | Signal handlers, parameter-path resolution, miscellaneous helpers |

### Calibration (Optional)

| File | Description |
|------|-------------|
| `calibration.c` / `calibration.h` | NLopt-driven automated calibration with multi-objective scoring |
| `calibration_stub.c` | No-op stub when calibration is disabled at compile time |
| `calibration_helpers.c` / `calibration_helpers.h` | Pure scoring, target summary, and timebase alignment utilities (always compiled) |
| `nlopt_winshim.c` | Windows compatibility shim for NLopt library |

## Key Design Principles

1. **Configuration-driven** — All calibratable parameters read from `params.txt`; no hardcoded physics overrides in source code.
2. **Literature-based formulations** — Every process equation cites peer-reviewed literature (Savenije, Fennel, Arndt, DOE, etc.) in code comments.
3. **Budget-preserving numerics** — Patankar (1980) implicit scheme for positivity; oxygen-availability cap prevents numerical undershoot; ordered dependency chain for metabolic safeguards (TOC → O₂ → NH₄ → NO₃ → PO₄ → Si).
4. **Transparent tropical adaptations** — Discharge-driven turbidity (`KD_flow`), salinity-stress phytoplankton mortality, quadratic closure for implicit zooplankton grazing.
5. **Strict validation** — Non-finite or significantly negative tracers trigger hard abort (`diagnostics_validate_tracer_state`), not silent clamping.

## Environment Variables (Runtime Overrides)

These are convenience switches for development/testing and do **not** alter model physics when unset:

| Variable | Purpose |
|----------|---------|
| `CGEM_TRANSPORT_ONLY` | Run hydrodynamics + transport without biogeochemistry |
| `CGEM_WARMUP_DAYS` | Override warmup period duration |
| `CGEM_MAXT_DAYS` | Override total simulation length |
| `CGEM_MASS_DIAG_DAYS` | Limit run to N days for quick mass-balance checks |
| `CGEM_FORCE_NO_CALIBRATION` | Disable calibration even if compiled with NLopt |
| `CGEM_SKIP_VALIDATION` | Skip post-simulation validation scripts |
| `CGEM_DEBUG_LEVEL` | Set diagnostic verbosity (0–3) |

For **production runs** (manuscript results), none of these variables are set — the model runs from `params.txt` configuration only.
