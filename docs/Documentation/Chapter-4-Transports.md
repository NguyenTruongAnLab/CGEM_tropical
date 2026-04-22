# Transport Module

C-GEM solves the 1D advection-dispersion equation (ADE) for each of the 16 chemical species. Transport is split into advection (TVD scheme) and dispersion (Crank-Nicolson), applied sequentially at each timestep.

## Governing Equation

$$\frac{\partial C}{\partial t} + U \frac{\partial C}{\partial x} = \frac{\partial}{\partial x}\left(D \frac{\partial C}{\partial x}\right) + S + R$$

| Symbol | Meaning | Unit |
|--------|---------|------|
| $C$ | Concentration | mg/L |
| $U$ | Cross-sectional velocity (from hydrodynamics) | m/s |
| $D$ | Longitudinal dispersion coefficient | m²/s |
| $S$ | External sources (tributaries, point sources) | mg/L/s |
| $R$ | Biogeochemical reaction rates | mg/L/s |

## Dispersion: Van der Burgh Framework

C-GEM parameterises dispersion using the Van der Burgh relationship (Gisen et al., 2015):

$$\frac{dD}{dx} = -K \cdot U \cdot h$$

where $K$ is the Van der Burgh coefficient (dimensionless, typically 0.2–0.8). Integrating for an exponentially converging estuary:

$$D(x) = D_0 \exp\!\left(-K \frac{x}{L_b}\right)$$

| Parameter | Meaning | Typical range |
|-----------|---------|---------------|
| $D_0$ | Dispersion at the mouth | 100–1000 m²/s |
| $K$ | Van der Burgh coefficient | 0.2–0.8 |
| $L_b$ | Width convergence length | site-specific |

These are calibrated in **Stage 1** alongside the Chézy coefficient, using salinity as the calibration target.

## Advection: TVD Scheme

The Total Variation Diminishing (TVD) scheme preserves sharp concentration fronts without introducing spurious oscillations.

**Algorithm (in `TVD()`, `transport.c`):**

1. Compute local CFL number from velocity and $\Delta x / \Delta t$
2. If CFL > 1, sub-step (up to 20 sub-steps per transport step)
3. For each cell interface, compute the upwind flux
4. Apply **Superbee** flux limiter to blend first-order and second-order fluxes
5. Update concentrations conservatively:

$$C_i^{n+1} = C_i^n - \frac{\Delta t}{A_i \Delta x}\left(F_{i+1/2} - F_{i-1/2}\right)$$

The Superbee limiter is defined as:

$$\phi(r) = \max\!\bigl(0,\;\min(2r,1),\;\min(r,2)\bigr)$$

where $r$ is the ratio of consecutive concentration gradients. This gives second-order accuracy in smooth regions and reduces to first-order upwind at sharp fronts.

## Dispersion: Crank-Nicolson

Dispersion is solved implicitly to avoid the strict $\Delta t < \Delta x^2 / (2D)$ stability limit that explicit treatment would impose:

$$\frac{C_i^{n+1} - C_i^n}{\Delta t} = \frac{1}{2}\left[\mathcal{L}_D(C^n) + \mathcal{L}_D(C^{n+1})\right]$$

where $\mathcal{L}_D$ is the dispersion operator. This forms a tridiagonal system solved by the Thomas algorithm (same solver used for hydrodynamics).

## Tributary and Source Inputs

Tributary inputs are injected as mass fluxes at specified cell indices:

$$\Delta C_i = \frac{Q_{\text{trib}} \cdot C_{\text{trib}}}{A_i \cdot \Delta x} \cdot \Delta t$$

Multiple tributaries can be configured in `config_input.txt`. Each tributary specifies:

- Cell index (injection location)
- Discharge time series
- Concentration time series for each species

## Boundary Conditions

| Boundary | Treatment |
|----------|-----------|
| Downstream (sea) | Prescribed concentration during flood; zero-gradient during ebb |
| Upstream (river) | Prescribed concentration at all times |

The flood/ebb switch uses the sign of the boundary velocity: when flow enters the domain, concentrations are imposed; when flow exits, the model uses the interior concentration (zero-gradient).

## Code Structure

| Function | File | Purpose |
|----------|------|---------|
| `Transport(t)` | `transport.c` | Main entry: loops over species |
| `TVD(co, s)` | `transport.c` | TVD advective fluxes |
| `Dispersion(co, s)` | `transport.c` | Crank-Nicolson dispersion |
| `bgboundary(t)` | `bcforcing.c` | Interpolate boundary concentrations |

Flux diagnostics (advective and dispersive flux per species per cell) are output to `OUT/Flux_Advection_{species}.csv` and `OUT/Flux_Dispersion_{species}.csv`.

## References

- Gisen, J.I.A., et al. (2015). *Predictive dispersion in alluvial estuaries.* HESS, 19, 2495-2514.
- Yee, H.C. (1987). *TVD schemes for nonlinear conservation laws.* NASA TM 89464.
