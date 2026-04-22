# Numerical Methods

C-GEM solves three coupled systems at each timestep — hydrodynamics, transport, and biogeochemistry — each with a numerical method suited to its physics.

## Solution Sequence

```
for each timestep (Δt = DELTI; bundled example = 360 s):
  1. Hyd(t)        — Saint-Venant equations (semi-implicit)
  2. bgboundary(t) — Interpolate boundary conditions
  3. Transport(t)  — Advection-dispersion for all species (TVD)
  4. Biogeo(t)     — Biogeochemical reactions (modified Euler)
  5. updateSPM(t)  — Sediment erosion/deposition
  6. Output(t)     — Write results (at configured intervals)
```

## Spatial Discretization

The domain is divided into $M$ cells with uniform spacing $\Delta x$ (default 2000 m). Variables are placed on a **staggered grid**:

| Position | Variables | Physical role |
|----------|-----------|---------------|
| Odd cells (1, 3, 5, ...) | Water level $\zeta$, concentrations $C_i$ | Scalar quantities |
| Even cells (2, 4, 6, ...) | Velocity $U$, discharge $Q$ | Flow quantities |

This staggering prevents checkerboard oscillations and naturally centres the pressure gradient at velocity nodes.

## Hydrodynamics: Semi-Implicit Scheme

Water levels are solved **implicitly** via a tridiagonal system; velocities are updated **semi-implicitly** using the new water levels. Friction is linearised using the previous-timestep velocity magnitude.

**Discretized continuity:**

$$\frac{h_i^{n+1} - h_i^n}{\Delta t} + \frac{Q_{i+1}^{n+1} - Q_{i-1}^{n+1}}{2\Delta x B_i} = 0$$

**Discretized momentum (solved for $U_i^{n+1}$):**

$$U_i^{n+1} = \frac{U_i^n - \Delta t\left[ U_i^n \frac{U_{i+1}^n - U_{i-1}^n}{2\Delta x} + g\frac{h_{i+1}^{n+1} - h_{i-1}^{n+1}}{2\Delta x} \right]}{1 + \Delta t\, \frac{g|U_i^n|}{C_i^2\, h_i^n}}$$

The resulting tridiagonal system $a_i h_{i-1} + b_i h_i + c_i h_{i+1} = d_i$ is solved by the Thomas algorithm in `NewUH()`. Unlike explicit schemes, this is unconditionally stable — no CFL restriction, allowing multi-minute timesteps; the bundled example currently uses `DELTI = 360 s`.

### Code path

`Hyd(t)` → `Newbc(t)` (boundaries) → `Coeffa(t)` (assemble matrix) → `NewUH(t)` / `Tridag()` (solve) → `Update()` (velocities) → `ApplyResults()`

## Transport: TVD Advection + Crank-Nicolson Dispersion

### Advection

The TVD scheme preserves sharp concentration fronts without introducing oscillations:

1. Compute CFL number; if CFL > 1, sub-step (up to 20 sub-steps)
2. For each cell, compute upwind flux and apply the **Superbee** flux limiter
3. Accumulate fluxes across sub-steps
4. Update concentrations conservatively: $C_i^{n+1} = C_i^n - \frac{\Delta t}{A_i \Delta x}(F_{i+1/2} - F_{i-1/2})$

### Dispersion

Dispersion is solved with a Crank-Nicolson (implicit) discretization, forming another tridiagonal system solved by the Thomas algorithm. This avoids the stability limit that explicit dispersion would impose.

### Code path

`Transport(t)` → for each species: `TVD(co, s)` (advective fluxes) → `Dispersion(co, s)` (implicit diffusion) → apply tributary sources

## Biogeochemistry: Modified Euler with Sub-stepping

Biogeochemical reactions are stiff (fast nitrification can have timescales of seconds). C-GEM uses a modified Euler (predictor-corrector) method with adaptive sub-stepping:

1. Compute reaction rates $R_i(C^n)$
2. Predictor: $C^* = C^n + \Delta t_{\text{bio}} R_i(C^n)$
3. Compute rates at predicted state: $R_i(C^*)$
4. Corrector: $C^{n+1} = C^n + \frac{\Delta t_{\text{bio}}}{2}[R_i(C^n) + R_i(C^*)]$
5. Repeat until the full `DELTI` transport step is covered

### Code path

`Biogeo(t)` → for each cell: `computePrimaryProduction()` → `computeBiogeochemicalReactions()` → `computeCarbonateChemistry()` → `updateBiogeochemicalState()`

## Mass Conservation

Mass conservation is enforced at every level:

| Component | Mechanism |
|-----------|-----------|
| Hydrodynamics | Implicit continuity equation with staggered grid |
| Advection | Conservative TVD flux form |
| Dispersion | Symmetric Crank-Nicolson stencil |
| Reactions | Stoichiometric coupling (Redfield ratios, elemental budgets) |
| Tributaries | Direct mass injection proportional to discharge |

The mass diagnostic (`mass_diag.c`) integrates total mass across all cells for each species and writes `OUT/mass_diag.csv` at configurable intervals, enabling verification of numerical conservation.

## Stability Diagnostics

`check_numerical_stability()` in `diagnostics.c` is called every timestep and aborts if it detects:

- NaN or Inf in any velocity, water level, or concentration
- Concentrations outside physically realistic bounds
- Velocity exceeding `MAX_VELOCITY`
