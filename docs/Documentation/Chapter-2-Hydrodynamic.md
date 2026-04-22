# Hydrodynamic Module

The hydrodynamics module solves the 1D Saint-Venant equations on a staggered grid using a semi-implicit predictor-corrector scheme. It provides velocity, water level, depth, and cross-sectional area fields that drive transport and biogeochemistry.

## Governing Equations

**Continuity (mass conservation):**

$$\frac{\partial \zeta}{\partial t} + \frac{1}{B r_s}\frac{\partial (BhU)}{\partial x} = 0$$

**Momentum:**

$$\frac{\partial U}{\partial t} + U\frac{\partial U}{\partial x} + g\frac{\partial \zeta}{\partial x} + \frac{g|U|U}{C^2 h} = 0$$

| Symbol | Meaning | Unit |
|--------|---------|------|
| $\zeta$ | Water surface elevation | m |
| $B$ | Channel width | m |
| $h$ | Water depth | m |
| $U$ | Cross-sectional average velocity | m/s |
| $C$ | Chézy coefficient | m¹ᐟ²/s |
| $r_s$ | Storage width ratio | — |
| $g$ | Gravitational acceleration (9.81) | m/s² |

### Momentum Terms

Each term in the momentum equation has a clear physical meaning:

| Term | Force | Dominance |
|------|-------|-----------|
| $\partial U/\partial t$ | Local acceleration (inertia) | During tidal reversals |
| $U \partial U/\partial x$ | Convective acceleration | In converging channels |
| $g \partial \zeta/\partial x$ | Pressure gradient | Primary driver of flow |
| $g\|U\|U/(C^2 h)$ | Bed friction | Always opposes flow |

Chézy is used rather than Manning because it is depth-independent and consistent with Savenije's estuarine theory.

## Estuarine Geometry

C-GEM uses exponential width convergence (Savenije, 2005):

$$B(x) = B_0 \exp(-x/L_b)$$

The estuary is divided into up to 4 segments, each with its own $B_0$, $L_b$, Chézy, and storage ratio. Depth is read from `INPUT/Geometry/river_bed_profile.txt` or computed from the analytical geometry. Width, depth, and cross-section area are computed at initialization and written to `OUT/geometry.csv`.

## Numerical Scheme

### Staggered Grid

The domain is discretized into $M$ cells with spacing $\Delta x$:

- **Odd cells (1, 3, 5, ...):** water levels $\zeta$
- **Even cells (2, 4, 6, ...):** velocities $U$

This staggering prevents checkerboard oscillations and ensures mass conservation.

### Semi-Implicit Time Integration

The discretized equations are:

**Continuity:**
$$\frac{h_i^{n+1} - h_i^n}{\Delta t} + \frac{Q_{i+1}^{n+1} - Q_{i-1}^{n+1}}{2\Delta x} = 0$$

**Momentum:**
$$\frac{Q_i^{n+1} - Q_i^n}{\Delta t} + \left[\frac{Q^2/A}\right]_{\text{explicit}} + gA_i^n \frac{h_{i+1}^{n+1} - h_{i-1}^{n+1}}{2\Delta x} + \frac{gQ_i^{n+1}|Q_i^n|}{C_i^2A_i^nR_i^n} = 0$$

Key: pressure gradient and friction are implicit → unconditionally stable. Advection is explicit. This allows large timesteps (3 min) without CFL restriction.

### Tridiagonal Solver

The implicit equations form a tridiagonal system:

$$a_i h_{i-1}^{n+1} + b_i h_i^{n+1} + c_i h_{i+1}^{n+1} = d_i$$

Solved by the Thomas algorithm (forward elimination + back-substitution).

### Code Implementation

| Step | Function | File |
|------|----------|------|
| Set boundary conditions | `Newbc(t)` | `hydrodynamics.c` |
| Assemble tridiagonal matrix | `Coeffa(t)` | `hydrodynamics.c` |
| Solve for water levels | `NewUH(t)` / `Tridag()` | `hydrodynamics.c` |
| Update velocities | `Update()` | `hydrodynamics.c` |
| Apply results | `ApplyResults()` | `hydrodynamics.c` |
| All steps combined | `Hyd(t)` | `hydrodynamics.c` |

### Boundary Conditions

| Boundary | Variable | Source |
|----------|----------|--------|
| Downstream (mouth) | Water level | `Tide(t)` — interpolated from observed elevation time series |
| Upstream (river) | Discharge | `Discharge_ups(t)` — interpolated from observed discharge |

## Tidal Diagnostics

C-GEM tracks tidal characteristics at each cell:

| Array | What it stores |
|-------|---------------|
| `velocity_min[]`, `velocity_max[]` | Velocity extremes per cell |
| `daily_min_level[]`, `daily_max_level[]` | Water level extremes |
| `tidal_range_mean[]` | Time-averaged tidal range |
| `dispersion_min[]`, `dispersion_max[]` | Dispersion bounds |

These are output to `OUT/` and can be used to diagnose tidal amplification, damping, and flood/ebb asymmetry.

## Stability Checks

The function `check_numerical_stability()` in `diagnostics.c` is called every timestep. It detects:

- NaN/Inf in velocity or water level
- Velocity exceeding `MAX_VELOCITY` (15 m/s)
- Depth below `MIN_DEPTH` (0.01 m) or above `MAX_DEPTH` (100 m)

If triggered, the model aborts with a diagnostic message.

## Calibration Parameters

Hydrodynamics are calibrated in **Stage 1** of the 4-stage workflow:

| Parameter | Typical range | Effect |
|-----------|--------------|--------|
| Chézy per segment | 30–70 m¹ᐟ²/s | Friction → tidal damping |
| Storage ratio $r_s$ | 1.0–1.5 | Tidal flat storage |
| `AMPL` | Site-specific | Tidal amplitude at mouth |

Lower Chézy → more friction → faster tidal damping upstream.

## References

- Savenije, H.H.G. (2005). *Salinity and Tides in Alluvial Estuaries.* Elsevier.
- Savenije, H.H.G. (2012). *Salinity and Tides in Alluvial Estuaries.* 2nd ed.
