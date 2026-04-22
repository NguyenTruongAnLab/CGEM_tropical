# Supplement 4 — Tributaries, Dispersion, and Diagnostics (Peer‑review note)

This document formalizes the recent transport improvements and diagnostics so they can be evaluated on scientific grounds and reproduced across estuaries with or without tributaries.

## 1. Scientific foundations

- Confluence mixing and tributary mass addition follow Fischer et al. (1979), ch. 7. We model tributaries as direct mass sources to the advection–dispersion equation; mixing and redistribution occur through the main solver. This prevents non-physical overwrites and guarantees accumulation toward a transport-driven steady state.
- Longitudinal dispersion follows the predictive framework of Gisen et al. (2015), anchored at the mouth by a physics-based estimate and propagated upstream by the Van der Burgh relationship (Savenije, 2012). This maintains a geometry/hydrodynamics-first philosophy and minimizes tunable parameters.

## 2. Tributary implementation (Fischer 1979)

We implement tributary inflows as mass source terms in the concentration equation for species C:

$$
\frac{\partial C}{\partial t} + \frac{\partial (uC)}{\partial x} = \frac{\partial}{\partial x}\left( D\, \frac{\partial C}{\partial x} \right) + S_C(x,t),\quad\text{with}\quad S_C=\frac{Q_{\mathrm{trib}}}{V_{\mathrm{eff}}}\big(C_{\mathrm{trib}}-C\big)\ \text{in the confluence cell.}
$$

- Per time step of duration $\Delta t$, the concentration increment added in the confluence control volume is

$$
\Delta C = \left(\frac{Q_{\mathrm{trib}}}{V_{\mathrm{eff}}}\right)\big(C_{\mathrm{trib}}-C\big)\,\Delta t.
$$

- Near‑field mixing volume. Following Fischer (1979, ch. 7) the near‑field confluence mixing zone has a longitudinal extent of a few channel widths. We represent this subgrid effect by an effective junction volume

$$
V_{\mathrm{eff}} = A\,L_{\mathrm{mix}}, \qquad L_{\mathrm{mix}} = \min\big(\Delta x,\ k_{\mathrm{mix}}\,B\big),\ \ k_{\mathrm{mix}}\in [2,5].
$$

Here $A$ is cross‑sectional area, $B$ width, and $\Delta x$ the grid spacing. The default $k_{\mathrm{mix}}=3$ can be considered a literature‑based constant; it should not be case‑tuned unless clear site evidence exists.

- Rationale. This preserves volume in the 1‑D main channel while correctly adding mass. The sign of the source is physically consistent: if $C_{\mathrm{trib}}<C$, then $S_C<0$ (dilution), and if $C_{\mathrm{trib}}>C$, then $S_C>0$ (enrichment). It prevents non‑physical “instant overwrite to $C_{\mathrm{trib}}$” and yields bounded accumulation across steps, approaching a steady transport balance near urban inputs.

- Non‑negativity. After transport we enforce $C\ge 0$ for all biogeochemical tracers; any negative value is clamped to zero to respect physical bounds.

Implementation
- File `src/transport.c`: `applyTributarySourceTerms()` computes $\Delta C$ and adds it to the confluence cell using $V_{\mathrm{eff}}$.
- Config (read from `INPUT/params.txt`): `enable_subgrid_confluence_mixing` (default 1), `junction_kmix_widths` (default 3.0 widths). Both are literature‑driven, not calibration targets.

## 3. Predictive dispersion (Gisen 2015 + Van der Burgh)

Mouth dispersion $D_0$ is predicted from hydrodynamics and geometry (Gisen et al., 2015, Eq. 41):

$$
\frac{D_0}{v\,E}=0.3958\,N_r^{0.57}\left(\frac{g}{C^2}\right)^{0.21},\qquad N_r=\frac{\Delta\rho/\rho\ \ g\ H\ Q_r\ T}{v^2\ A\ E}.
$$

- $v$ is the tidal velocity amplitude (we use the first velocity interface), $E = v\,T/\pi$ is the tidal excursion, $T$ is the M2 period, $C$ is the Chezy coefficient, $g$ gravity, $H$ depth, $A$ cross‑sectional area, $Q_r$ the net river discharge magnitude, and $\Delta\rho/\rho$ stratification parameter (Savenije, 2012). Scalar checks and floors avoid division by zero.
- Upstream profile. We map $D_0$ inland using the Van der Burgh relationship tied to geometry:

$$
D(x)=D_0\left( \frac{A(x)}{A_0} \right)^{K_{\mathrm{VDB}}},\qquad K_{\mathrm{VDB}}\in[0.5,0.8].
$$

- Safety rails. We constrain $D\in[\text{MIN\_DISP\_COEFF},\text{MAX\_DISP\_COEFF}]$ to prevent unphysical extremes; these are not calibration targets but numerical guards.

Implementation
- File `src/transport.c`: `Dispcoef()` implements the equations above. Parameters are configuration‑driven: `K_VDB_ENHANCED` (the single recommended tuning knob), `DELTA_RHO_OVER_RHO_PARAM` (fixed 0.02–0.03), and optional `MIN/MAX` clamps.

Minimal tunable set
- Predictive‑first philosophy: tune at most one number, $K_{\mathrm{VDB}}$. The density ratio and clamps should remain fixed. An optional global multiplier (`enhancement_factor`) can be set to 1.0 to remove it entirely.

## 4. Warmup and boundary forcing consistency

- Time series policy. Hourly series (elevation, light) remain dynamic during warmup; daily series (upstream net discharge, tributary discharges and concentrations, chemical BCs) are held constant at “day 1” for steady transport equilibration. Thus, tributaries inject mass continuously during warmup with stable daily values—no shock at the warmup/post‑warmup boundary.
- Auto‑equilibration. Optionally extend warmup until transport reaches steady state: compute relative day‑to‑day change of domain‑mean (NH4, NO3, Sal) and require $\max\Delta\le \varepsilon$ (default $1\%$) for $N$ consecutive days (default 3). This prevents starting reactive simulations from a non‑equilibrated transport state.

## 5. Mass‑balance diagnostics with tributary sources

Let total domain mass for species $C$ be $M=\sum_i C_i V_i$. Over a validation interval $\Delta t_v$ we estimate

$$
\frac{\mathrm{d}M}{\mathrm{d}t}\frac{M(t)-M(t-\Delta t_v)}{\Delta t_v}.
$$

We now account for internal source terms from tributaries:

$$
R = \frac{\mathrm{d}M}{\mathrm{d}t} - S_{\mathrm{trib}},\qquad S_{\mathrm{trib}}=\sum_k Q_{\mathrm{trib},k}\,C_{\mathrm{trib},k},
$$

with $R$ the residual. In reactions‑off mode, $R$ should be small when boundary conditions are consistent. The diagnostic reports $\mathrm{d}M/\mathrm{d}t$, $S_{\mathrm{trib}}$, and the residual.

Implementation
- File `src/transport.c`: accumulates per‑species $S_{\mathrm{trib}}$ each step.
- File `src/diagnostics.c`: `validate_mass_conservation()` subtracts $S_{\mathrm{trib}}$ for the tested species (e.g., NH4). Salinity checks remain unaffected (no internal source).

## 6. Configuration summary (minimal, literature‑based)

- Subgrid confluence mixing (Fischer 1979):
  - `enable_subgrid_confluence_mixing` (default 1)
  - `junction_kmix_widths` (default 3.0 widths; literature 2–5; treat as fixed unless field evidence suggests otherwise)
- Predictive dispersion (Gisen 2015 + Van der Burgh):
  - `K_VDB_ENHANCED` (0.5–0.8): recommended single tuning knob
  - `DELTA_RHO_OVER_RHO_PARAM` (0.02–0.03): fixed physics
  - `MIN_DISP_COEFF`, `MAX_DISP_COEFF`: safety rails (not calibration)
  - `enhancement_factor`: set 1.0 to remove; optional otherwise

## 7. Edge cases and guarantees

- Non‑negativity: $C\ge 0$ enforced after each transport step.
- No tributaries: the source term is simply zero; transport remains predictive.
- Stability: junction source is mass‑only; hydrodynamic volume is preserved; $L_{\mathrm{mix}}\le \Delta x$ prevents nonlocal updates.
- Time step: no $\Delta t$ change is required.

## 8. Validation protocol

- Phase‑1 (reactions OFF): verify accumulation at urban confluences in steady warmup, then persistence along the main stem. NH4 budget: check $\mathrm{d}M/\mathrm{d}t$, $S_{\mathrm{trib}}$, residual.
- Sensitivity (optional): if peaks are overly smoothed, adjust only $K_{\mathrm{VDB}}$ within 0.5–0.8; keep other parameters fixed.

## 9. References

- Fischer, H. B., List, E. J., Koh, R. C. Y., Imberger, J., & Brooks, N. H. (1979). Mixing in Inland and Coastal Waters. Academic Press. (Confluences: Chapter 7)
- Savenije, H. H. G. (2012). Salinity and Tides in Alluvial Estuaries. Completely Revised 2nd Edition. Delft Academic Press.
- Gisen, J. I. A., Savenije, H. H. G., & Nijzink, R. C. (2015). Revised predictive equations for salt intrusion modelling in estuaries. Hydrol. Earth Syst. Sci., 19, 2791–2803.