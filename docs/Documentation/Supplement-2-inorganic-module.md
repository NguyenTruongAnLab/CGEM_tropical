# Supplement 2 — Inorganic Carbon Module

This supplement documents the carbonate chemistry implementation in C-GEM, including equilibrium constants, the pH solver, and CO₂ gas exchange. The module computes pCO₂, pH, and CO₂ diagnostically from transported DIC and AT at each timestep.

## Carbonate Equilibria

$$\text{CO}_2(\text{aq}) \xrightleftharpoons{K_1} \text{H}^+ + \text{HCO}_3^- \xrightleftharpoons{K_2} 2\text{H}^+ + \text{CO}_3^{2-}$$

Given DIC and [H⁺], the species fractions are:

$$[\text{CO}_2] = \frac{\text{DIC}}{1 + K_1/[\text{H}^+] + K_1 K_2/[\text{H}^+]^2}$$

$$[\text{HCO}_3^-] = \frac{\text{DIC}}{1 + [\text{H}^+]/K_1 + K_2/[\text{H}^+]}$$

$$[\text{CO}_3^{2-}] = \frac{\text{DIC}}{1 + [\text{H}^+]/K_2 + [\text{H}^+]^2/(K_1 K_2)}$$

## Equilibrium Constants

All constants are temperature ($T$ in K) and salinity ($S$) dependent.

| Constant | Function in code | Reference |
|----------|-----------------|-----------|
| $K_1$ (1st dissoc.) | `disK1(Sal, t)` | Cai & Wang (1998) |
| $K_2$ (2nd dissoc.) | `disK2(Sal, t)` | Cai & Wang (1998) |
| $K_w$ (water) | `disK3(Sal, t)` | Dickson & Riley (1979) |
| $K_B$ (boric acid) | `disK5(Sal, t)` | Dickson (1990) |
| $K_0$ (CO₂ solubility) | `calc_co2_solubility(Tk, S)` | Weiss (1974) |

**$K_1$ formulation:**

```c
double pK1 = -14.8425 + (3404.71/T) + (0.032786*T);
double f1 = -0.0230848 - (14.3456/T);
double f2 = 0.000691881 + (0.429955/T);
K1 = pow(10, -(pK1 + f1*sqrt(Sal) + f2*Sal));
```

**$K_0$ (Henry's Law) formulation (Weiss, 1974):**

```c
double i0 = -58.0931 + (90.5069/Tk) + 22.2940*log(Tk);
i0 += Sal * (0.027766 - 0.025888*(Tk/100.0) + 0.0050578*pow(Tk/100.0, 2));
K0 = exp(i0);  // mol kg⁻¹ atm⁻¹
```

## pH Solver

The solver finds [H⁺] iteratively from known DIC and AT:

1. Isolate carbonate alkalinity: $\text{AT}_\text{carb} = \text{AT} - [\text{B(OH)}_4^-] - [\text{OH}^-] + [\text{H}^+]$
2. Semi-analytical update:

$$[\text{H}^+]_\text{new} = \frac{1}{2}\left[\left(\frac{\text{DIC}}{\text{AT}_\text{carb}} - 1\right) K_1 + \sqrt{\left(1 - \frac{\text{DIC}}{\text{AT}_\text{carb}}\right)^2 K_1^2 - 4K_1 K_2\left(1 - \frac{2\,\text{DIC}}{\text{AT}_\text{carb}}\right)}\right]$$

3. Recalculate species fractions and repeat until $|\text{AT}_\text{calc} - \text{AT}| < \text{TOL}$ (max 50 iterations)
4. Convert: $\text{pH} = -\log_{10}[\text{H}^+]$

Total boron is computed from salinity: $t_b = 416.0 \times (S / 35)$ µmol/kg.

## Biogeochemical Effects on DIC and Alkalinity

| Process | $\Delta$DIC / mol C | $\Delta$AT / mol C | Sign convention |
|---------|---------------------|---------------------|-----------------|
| NPP via NO₃ | −1 | +17/106 | Uptake = negative DIC |
| NPP via NH₄ | −1 | −15/106 | |
| Aerobic respiration | +1 | +15/106 | Release = positive DIC |
| Nitrification | 0 | −2 per mol N | Carbon-neutral |
| Denitrification | +1 | +93.4/106 | |

Code implementation:

```c
reactionDIC[i] = adegrad[i] + denit[i] - NPP[i];

reactionTA[i] = (15.0/106.0) * adegrad[i]
              + (93.4/106.0) * denit[i]
              - 2.0 * nitrif[i]
              + (-15.0/106.0) * (NPP_NH4_Phy1 + NPP_NH4_Phy2)
              + (17.0/106.0) * (NPP_NO3_Phy1 + NPP_NO3_Phy2);
```

Stoichiometric coefficients follow Wolf-Gladrow et al. (2007) and Soetaert et al. (2007).

## CO₂ Gas Exchange

$$F_{\text{CO}_2} = \frac{k_{\text{CO}_2}}{h} \cdot \bigl([\text{CO}_2] - K_0 \cdot x_{\text{CO}_2} \cdot \rho_w\bigr)$$

The gas transfer velocity combines flow and wind components:

| Component | Formula | Reference |
|-----------|---------|-----------|
| Flow-induced | $k_\text{flow} = \sqrt{U \cdot D_m / h}$ | O'Connor & Dobbins (1958) |
| Wind-induced | $k_\text{wind} = 0.31\, U_{10}^2 \cdot (Sc/660)^{-0.5}$ | Wanninkhof (1992) |

Schmidt number for CO₂:

$$Sc = 2073.1 - 125.62\,T + 3.6276\,T^2 - 0.043219\,T^3$$

The total velocity is scaled by calibration parameter `piston_velocity_scale` (Stage 4).

## Processing Sequence

For each cell at each timestep:

1. Compute biogeochemical $\Delta$DIC and $\Delta$AT
2. Update DIC and AT
3. Solve pH iteratively (warm-started from previous timestep)
4. Compute CO₂, HCO₃⁻, CO₃²⁻ from DIC + pH
5. Compute CO₂ saturation from $K_0$, temperature, salinity
6. Compute gas transfer velocity $k$
7. Compute CO₂ flux and update DIC
8. Compute pCO₂ = [CO₂] / $K_H$

## References

- Cai, W.-J. & Wang, Y. (1998). *Limnol. Oceanogr.*, 43, 657-668.
- Dickson, A.G. & Riley, J.P. (1979). *Marine Chemistry*, 7, 89-99.
- Dickson, A.G. (1990). *Deep-Sea Research*, 37, 755-766.
- Weiss, R.F. (1974). *Marine Chemistry*, 2, 203-215.
- Wanninkhof, R. (1992). *J. Geophys. Res.*, 97, 7373-7382.
- Wolf-Gladrow, D.A. et al. (2007). *Marine Chemistry*, 106, 287-300.
- Soetaert, K. et al. (2007). *Science of the Total Environment*, 373, 62-73.
