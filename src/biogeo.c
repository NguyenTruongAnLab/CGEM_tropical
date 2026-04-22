/**
 * @file biogeo.c
 * @brief Biogeochemical reaction module for C-GEM
 *
 * @details
 * This file implements the biogeochemical reaction network
 *
 * @author An Nguyen
 * @date 2025-11-25
 */

#include "define.h"
#include "variables.h"
#include "diagnostics.h"
#include "file.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "biogeo.h"
#include "utilities.h"

/* ========================================================================
 * CONSTANTS (Matching Fortran cgempars.f90)
 * ======================================================================== */
static const double EULER_CONSTANT = 0.5772156649;
static const double PH_TOL = 1e-10;
static const int PH_MAX_ITER = 50;

/**
 * @brief Efficient polynomial evaluation using Horner's method.
 *
 * Approximates the integral terms for small p (p <= 1.0).
 * Formula: p - p^2/4 + p^3/18 - p^4/96 + p^5/600
 */
static inline double poly_approx(double p) {
    return p * (1.0 - p * (0.25 - p * (1.0 / 18.0 - p * (1.0 / 96.0 - p / 600.0))));
}

/**
 * @brief Light limitation integral (Jassby & Platt 1976)
 * @note Direct port from Fortran CGEM_Functions.f90
 */
static double light_limitation(double alpha_phy, double Pbmax_val, double I0_val, 
                               double KD, double depth) {
    if (I0_val <= 0.0 || depth <= 0.0 || Pbmax_val <= 0.0) return 0.0;
    
    double Ebottom = I0_val * exp(-KD * depth);
    
    if (Ebottom < 1e-30) {
        // Asymptotic limit as Ebottom -> 0
        double psurf = I0_val * alpha_phy / Pbmax_val;
        if (psurf <= 0.0) return 0.0;
        if (psurf <= 1.0) {
            return poly_approx(psurf) / KD;
        } else {
            double E1_surf = exp(-psurf) / (psurf + 1.0 - 1.0 / (psurf + 3.0 - 
                             4.0 / (psurf + 5.0 - 9.0 / (psurf + 7.0 - 16.0 / (psurf + 9.0)))));
            return (E1_surf + EULER_CONSTANT + log(psurf)) / KD;
        }
    }
    
    double psurf = I0_val * alpha_phy / Pbmax_val;
    double pbot = Ebottom * alpha_phy / Pbmax_val;
    
    double appGAMMAsurf, appGAMMAbot;
    
    // Approximation of incomplete gamma function (exact Fortran formulas)
    if (psurf <= 1.0) {
        appGAMMAsurf = -(log(psurf) + EULER_CONSTANT - poly_approx(psurf));
    } else {
        appGAMMAsurf = exp(-psurf) / (psurf + 1.0 - 1.0 / (psurf + 3.0 - 
                       4.0 / (psurf + 5.0 - 9.0 / (psurf + 7.0 - 16.0 / (psurf + 9.0)))));
    }
    
    if (pbot <= 1.0) {
        appGAMMAbot = -(log(pbot) + EULER_CONSTANT - poly_approx(pbot));
    } else {
        appGAMMAbot = exp(-pbot) / (pbot + 1.0 - 1.0 / (pbot + 3.0 - 
                      4.0 / (pbot + 5.0 - 9.0 / (pbot + 7.0 - 16.0 / (pbot + 9.0)))));
    }
    
    return (appGAMMAsurf - appGAMMAbot + log(I0_val / Ebottom)) / KD;
}

/* ========================================================================
 * PISTON VELOCITY (Matching Fortran PistonVelocity function)
 * ======================================================================== */

/**
 * @brief Calculate piston velocity for gas exchange
 * @note Combines flow-induced and wind-induced components
 */
static double piston_velocity(int i, int t) {
    double vel = fabs(velocity[i]);
    double depth = waterDepth[i];
    if (depth <= 0.1) return 0.0;
    
    // Flow-induced component (O'Connor & Dobbins 1958)
    double kflow = sqrt(vel * Diff(t) / depth);
    
    // Wind-induced component (Wanninkhof 1992)
    double ws_val = windspeed(t, i);
    double Sc_val = Sc(t, i);
    double kwind = (1.0 / 3.6e5) * 0.31 * ws_val * ws_val * pow(Sc_val / 660.0, -0.5);
    
    double k = kflow + kwind;
    if (piston_velocity_scale > 0.0 && isfinite(piston_velocity_scale)) {
        k *= piston_velocity_scale;
    }
    return k;
}

// The hardcoded string-matching is_urban_tributary() has been removed.
// We strictly use the configuration-driven Tributary->is_urban flag now.

// Forward declarations for external functions

/* ========================================================================
 * CARBONATE SYSTEM (Robust, unit-consistent)
 * ======================================================================== */

/**
 * @brief Compute total alkalinity (TA) from DIC and H+ using a simplified carbonate system.
 *
 * @details
 * Uses the standard alkalinity definition (DOE 1994 / Dickson et al. 2007):
 *
 * $$TA = [HCO_3^-] + 2[CO_3^{2-}] + [B(OH)_4^-] + [OH^-] - [H^+] + [NH_4^+] - [NO_3^-]$$
 *
 * with carbonate speciation given by $K_1, K_2$ and borate by $K_B$.
 *
 * Units: all concentrations and equilibrium constants are treated on a mol/L scale.
 * This is critical because the model state variables DIC/AT are stored as mmol/m^3
 * (numerically equivalent to µmol/L), while pH is defined using [H+] in mol/L.
 */
static double carbonate_ta_from_dic_h_molL(
    double DIC_molL,
    double H_molL,
    double NH4_molL,
    double NO3_molL,
    double K1,
    double K2,
    double Kw,
    double TotB_molL,
    double Kb
) {
    if (!(DIC_molL >= 0.0) || !(H_molL > 0.0)) return NAN;

    const double H = H_molL;
    const double denom = 1.0 + (K1 / H) + (K1 * K2) / (H * H);
    if (!(denom > 0.0)) return NAN;

    const double CO2 = DIC_molL / denom;
    const double HCO3 = (K1 * CO2) / H;
    const double CO3 = (K2 * HCO3) / H;
    const double BOH4 = (TotB_molL > 0.0) ? (TotB_molL * Kb) / (Kb + H) : 0.0;
    const double OH = Kw / H;

    // Nutrient alkalinity terms included for consistency with the model's TA state updates
    // (Soetaert et al. 2007 style TA bookkeeping): +NH4 - NO3.
    return HCO3 + 2.0 * CO3 + BOH4 + OH - H + NH4_molL - NO3_molL;
}

/**
 * @brief Solve carbonate speciation from DIC and TA.
 *
 * @param dic_mmol_m3 DIC in mmol/m^3
 * @param ta_mmol_m3  Total alkalinity (equivalents) in mmol/m^3
 * @param sal         Salinity (psu)
 * @param t           Time (seconds)
 * @param i           Cell index (for temperature/salinity-dependent Henry)
 * @param pH_out      Output pH (dimensionless)
 * @param co2_mmol_m3_out Output aqueous CO2* in mmol/m^3
 * @param pco2_uatm_out Output pCO2 in µatm
 * @return 1 if solved successfully, 0 otherwise
 */
static int carbonate_solve_dic_ta(
    double dic_mmol_m3,
    double ta_mmol_m3,
    double nh4_mmol_m3,
    double no3_mmol_m3,
    double sal,
    int t,
    int i,
    double *pH_out,
    double *co2_mmol_m3_out,
    double *pco2_uatm_out
) {
    if (!pH_out || !co2_mmol_m3_out || !pco2_uatm_out) return 0;
    if (!(dic_mmol_m3 > 0.0) || !(ta_mmol_m3 > 0.0) || !isfinite(dic_mmol_m3) || !isfinite(ta_mmol_m3)) {
        return 0;
    }

    // Convert model units (mmol/m^3 == µmol/L) -> mol/L for carbonate algebra.
    const double DIC = dic_mmol_m3 * 1e-6;
    const double TA = ta_mmol_m3 * 1e-6;
    const double NH4 = fmax(nh4_mmol_m3, 0.0) * 1e-6;
    const double NO3 = fmax(no3_mmol_m3, 0.0) * 1e-6;

    // Equilibrium constants (same literature formulations as biogeout.c).
    const double K1 = disK1(sal, t);
    const double K2 = disK2(sal, t);
    const double Kw = disK3(sal, t);
    const double Kb = disK5(sal, t);

    // Total borate from DOE94: 416 µmol/kg at S=35. Approximate as mol/L.
    const double TotB = 416e-6 * (sal / 35.0);

    // Bracket pH in a broad but physically plausible range.
    double pH_lo = 4.0;
    double pH_hi = 10.5;
    double H_lo = pow(10.0, -pH_lo);
    double H_hi = pow(10.0, -pH_hi);

    double f_lo = carbonate_ta_from_dic_h_molL(DIC, H_lo, NH4, NO3, K1, K2, Kw, TotB, Kb) - TA;
    double f_hi = carbonate_ta_from_dic_h_molL(DIC, H_hi, NH4, NO3, K1, K2, Kw, TotB, Kb) - TA;
    if (!isfinite(f_lo) || !isfinite(f_hi)) return 0;

    // If not bracketed, try expanding slightly (rare; keeps solver robust).
    if (f_lo * f_hi > 0.0) {
        pH_lo = 3.0;
        pH_hi = 11.0;
        H_lo = pow(10.0, -pH_lo);
        H_hi = pow(10.0, -pH_hi);
        f_lo = carbonate_ta_from_dic_h_molL(DIC, H_lo, NH4, NO3, K1, K2, Kw, TotB, Kb) - TA;
        f_hi = carbonate_ta_from_dic_h_molL(DIC, H_hi, NH4, NO3, K1, K2, Kw, TotB, Kb) - TA;
        if (!isfinite(f_lo) || !isfinite(f_hi) || f_lo * f_hi > 0.0) {
            return 0;
        }
    }

    // Bisection in pH-space (monotonic in practice for DIC/TA > 0).
    double pH_mid = 0.5 * (pH_lo + pH_hi);
    double f_mid = 0.0;
    for (int iter = 0; iter < 40; iter++) {
        pH_mid = 0.5 * (pH_lo + pH_hi);
        const double H_mid = pow(10.0, -pH_mid);
        f_mid = carbonate_ta_from_dic_h_molL(DIC, H_mid, NH4, NO3, K1, K2, Kw, TotB, Kb) - TA;
        if (!isfinite(f_mid)) return 0;

        // Convergence: absolute TA error on mol/L scale.
        if (fabs(f_mid) < 1e-12) break;

        if (f_lo * f_mid <= 0.0) {
            pH_hi = pH_mid;
            f_hi = f_mid;
        } else {
            pH_lo = pH_mid;
            f_lo = f_mid;
        }
    }

    const double H_final = pow(10.0, -pH_mid);
    const double denom = 1.0 + (K1 / H_final) + (K1 * K2) / (H_final * H_final);
    if (!(denom > 0.0)) return 0;
    const double CO2_molL = DIC / denom;
    if (!isfinite(CO2_molL) || CO2_molL < 0.0) return 0;

    const double Henry = KH(t, i); // mmol/m^3.atm
    const double CO2_mmol_m3 = CO2_molL * 1e6;
    double pco2 = NAN;
    if (Henry > 1e-12) {
        pco2 = (CO2_mmol_m3 / Henry) * 1e6; // µatm
    }
    if (!isfinite(pco2) || pco2 < 0.0) return 0;

    *pH_out = pH_mid;
    *co2_mmol_m3_out = CO2_mmol_m3;
    *pco2_uatm_out = pco2;
    return 1;
}

/* ========================================================================
 * SEDIMENT FLUXES (Matching Fortran CGEM_Sediment.f90)
 * ======================================================================== */

/**
 * @brief Compute sediment erosion and deposition fluxes
 *
 * Erosion: Partheniades (1965) — E = M₀(τ_b/τ_ce − 1) when τ_b > τ_ce
 * Deposition: Krone (1962) — D = w_s(1 − τ_b/τ_cd)C when τ_b < τ_cd
 * Bottom shear: τ_b = ρgU²/C² (Chézy-based)
 */
static void compute_sediment_fluxes(int i) {
    if (i < 1 || i > M) return;
    
    double depth = waterDepth[i];
    double chezy_val = Chezy[i];
    
    if (depth <= 0.0 || chezy_val <= 0.0) {
        tau_b[i] = erosion_s[i] = deposition_s[i] = 0.0;
        erosion_v[i] = deposition_v[i] = 0.0;
        return;
    }
    
    // Hydrodynamics uses a staggered grid: velocities are defined on even indices
    // (faces), while state variables (including SPM) live on cell/control-volume
    // indices that are often odd. Using velocity[i] directly for odd i can yield
    // near-zero shear stress, suppressing erosion and exaggerating deposition.
    double vel = 0.0;
    if ((i % 2) == 0) {
        vel = velocity[i];
    } else {
        // Cell-centered representative velocity from adjacent faces.
        double vl = 0.0;
        double vr = 0.0;

        if (i - 1 >= 2) {
            vl = velocity[i - 1];
        } else if (i + 1 <= M1) {
            vl = velocity[i + 1];
        }

        if (i + 1 <= M1) {
            vr = velocity[i + 1];
        } else if (i - 1 >= 2) {
            vr = velocity[i - 1];
        }

        if (!isfinite(vl)) vl = 0.0;
        if (!isfinite(vr)) vr = 0.0;
        vel = 0.5 * (vl + vr);
    }
    
    // Bottom shear stress (τ_b = ρ g U² / C²)
    tau_b[i] = rho_w * g * vel * vel / (chezy_val * chezy_val);
    
    // Erosion rate (Partheniades 1965)
    const double Mero_eff = Mero[i]; 

    if (tau_b[i] > tau_ero[i] && tau_ero[i] > 0.0) {
        erosion_s[i] = Mero_eff * (tau_b[i] / tau_ero[i] - 1.0);
    } else {
        erosion_s[i] = 0.0;
    }
    
    // Deposition rate (Krone 1962)
    if (tau_dep[i] >= tau_b[i] && tau_dep[i] > 0.0) {
        deposition_s[i] = ws * (1.0 - tau_b[i] / tau_dep[i]) * v[SPM].c[i];
    } else {
        deposition_s[i] = 0.0;
    }
    
    // Convert surface fluxes to volumetric rates
    erosion_v[i] = erosion_s[i] / depth;
    deposition_v[i] = deposition_s[i] / depth;
}

/**
 * @brief Update suspended sediment concentrations (standalone for fast calibration)
 *
 * Standard Partheniades (1965) erosion and Krone (1962) deposition with
 * semi-implicit (Patankar) time stepping to guarantee positivity.
 */
void updateSuspendedSediment(int t) {
    if (!enable_suspended_sediment_dynamics) return;

    if (M <= 0) return;
    
    const double dt = (double)DELTI;
    
    for (int i = 1; i <= M; i++) {
        compute_sediment_fluxes(i);
        
        double C_old = v[SPM].c[i];
        double erosion_rate = erosion_v[i];
        double deposition_flux_explicit = deposition_v[i];
        
        // Infer specific deposition rate K_dep [1/s]
        double K_dep = 0.0;
        if (C_old > 1e-20) {
            K_dep = deposition_flux_explicit / C_old;
        }
        
        // Patankar semi-implicit update
        double C_new = (C_old + erosion_rate * dt) / (1.0 + K_dep * dt);
        
        v[SPM].c[i] = C_new;

        if (!isfinite(v[SPM].c[i])) {
              diagnostics_fatal_tracer_error("updateSuspendedSediment", t, SPM, i, v[SPM].c[i], "nonfinite");
              return;
        }
        if (v[SPM].c[i] < 0.0) {
               diagnostics_fatal_tracer_error("updateSuspendedSediment", t, SPM, i, v[SPM].c[i], "negative");
               return;
        }
    }
    
    if (M >= 1) v[SPM].c[0] = v[SPM].c[1];
}

/* ========================================================================
 * MAIN BIOGEOCHEMISTRY FUNCTION (Clean refactored version)
 * ======================================================================== */

/* ================================================================
     * PHOSPHORUS ADSORPTION HELPER
     * ================================================================ */
static void compute_phosphorus_adsorption_step(int i, double dt) {
        double SPM_kg = v[SPM].c[i];  // SPM is stored in kg/m³
        
        // Thresholds to avoid singularities/numerical noise
        if (SPM_kg > 1e-6 && v[PO4].c[i] > 1e-9) {
            // Langmuir equilibrium PIP concentration
            // PIP_eq = Pac * SPM * (PO4 / (PO4 + Kps))
            double PIP_eq = P_ac * SPM_kg * (v[PO4].c[i] / (v[PO4].c[i] + K_ps));
            
            // Kinetic relaxation: dPIP/dt = k_ads * (PIP_eq - PIP)
            // Explicit: PIP_new = PIP + k * dt * (PIP_eq - PIP)  --> Unstable if k*dt > 1
            // Implicit (Robust): PIP_new = (PIP + k * dt * PIP_eq) / (1 + k * dt)
            // This is unconditionally stable and non-oscillatory.
            
            double k_dt = k_ads * dt;
            double PIP_old = v[PIP].c[i];
            
            // Calculate target PIP_new based on relaxation logic
            double PIP_target = (PIP_old + k_dt * PIP_eq) / (1.0 + k_dt);
            double dPIP = PIP_target - PIP_old;

            // Mass constraint: Cannot adsorb more PO4 than exists
            // (Desorption is always safe because PIP_target is weighted average of PIP_old and PIP_eq)
            if (dPIP > 0.0 && dPIP > v[PO4].c[i]) {
                 dPIP = v[PO4].c[i]; // Physical limit of available reagent
            }
            // Note: We do not need artificial 'fmin/fmax' clamps for stability anymore, 
            // only this physical mass availability check.

            v[PIP].c[i] += dPIP;
            v[PO4].c[i] -= dPIP;
        }
        
        // P burial
        // dPIP/dt = - (Deposition / SPM) * PIP
        // Explicit: PIP_new = PIP - Rate * PIP * dt --> Negative if Rate*dt > 1
        // Implicit: PIP_new = PIP / (1 + Rate * dt) --> Always Positive
        
        if (v[PIP].c[i] > 1e-10 && deposition_v[i] > 0.0 && SPM_kg > 1e-6) {
            double specific_burial_rate = deposition_v[i] / SPM_kg; // [1/s]
            double burial_factor = 1.0 / (1.0 + specific_burial_rate * dt);
            
            double PIP_before = v[PIP].c[i];
            v[PIP].c[i] = PIP_before * burial_factor;
            
            // Flux tracking (optional, but good for diagnostics)
            // double mass_buried = (PIP_before - v[PIP].c[i]); 
        }
    }

/**
 * @brief Main biogeochemistry calculation (clean version)
 * 
 * This function consolidates all biogeochemical calculations into a single
 * efficient loop, following the structure of the Fortran CGEM_Biogeochemistry
 * subroutine while maintaining C-GEM's scientific foundation.
 * 
 * @param t Current simulation time (seconds)
 */
void Biogeo(int t) {

    // Fail-fast on invalid tracer state before any reactions (enables unit death tests).
    diagnostics_validate_tracer_state("Biogeo:entry", t);

    // Carbonate diagnostics health counters (throttled to daily summaries)
    static long carbonate_day = -1;
    static int carbonate_solver_failures_day = 0;
    static int carbonate_nonfinite_day = 0;
    static int carbonate_extreme_day = 0;

    long this_day = (long)((double)t / 86400.0);
    if (carbonate_day < 0) {
        carbonate_day = this_day;
    }
    if (this_day != carbonate_day) {
        if (debug_level > DEBUG_LEVEL_ESSENTIAL && (carbonate_solver_failures_day > 0 || carbonate_nonfinite_day > 0 || carbonate_extreme_day > 0)) {
            printf("🧪 Carbonate diagnostics (day %ld): solver_fail=%d, nonfinite=%d, extreme=%d\n",
                   carbonate_day, carbonate_solver_failures_day, carbonate_nonfinite_day, carbonate_extreme_day);
        }
        carbonate_day = this_day;
        carbonate_solver_failures_day = 0;
        carbonate_nonfinite_day = 0;
        carbonate_extreme_day = 0;
    }
    // Transport-only mode: no reactions. Optionally still compute carbonate
    // speciation diagnostics (pH/CO2/pCO2) from transported DIC+AT so Phase3
    // validators and downstream diagnostics have consistent carbonate fields.
    // This does NOT modify DIC/AT (mass-conservative).
    if (!enable_biogeochemical_reactions) {
        if (enable_carbonate_diagnostics) {
            for (int i = 1; i <= M; i++) {
                double sal = v[Sal].c[i];
                double sal_safe = (sal < 1e-10) ? 0.0 : sal;
                double pH_now = v[PH].c[i];
                double CO2aq_eq = v[CO2].c[i];
                double pCO2_now = v[pCO2].c[i];
                int carbonate_ok = carbonate_solve_dic_ta(v[DIC].c[i], v[AT].c[i], v[NH4].c[i], v[NO3].c[i], sal_safe, t, i, &pH_now, &CO2aq_eq, &pCO2_now);
                if (!carbonate_ok) {
                    carbonate_solver_failures_day++;
                    if (carbonate_solver_failure_policy == CARBONATE_FAIL_FATAL) {
                        fprintf(stderr,
                                "FATAL: Carbonate solver failed (transport-only diagnostics) at cell=%d t=%d (DIC=%.17g AT=%.17g Sal=%.17g)\n",
                                i, t, v[DIC].c[i], v[AT].c[i], sal_safe);
                        exit(EXIT_FAILURE);
                    }
                    continue;
                }
                if (!isfinite(pH_now) || !isfinite(CO2aq_eq) || !isfinite(pCO2_now)) {
                    carbonate_nonfinite_day++;
                    continue;
                }
                // Flag extreme pCO2 as a diagnostic signal (but do not alter mass fields)
                if (pCO2_now > 1.0e6) {
                    carbonate_extreme_day++;
                }
                v[PH].c[i] = pH_now;
                v[CO2].c[i] = CO2aq_eq;
                v[pCO2].c[i] = pCO2_now;
            }
            if (M >= 1) {
                v[PH].c[0] = v[PH].c[1];
                v[CO2].c[0] = v[CO2].c[1];
                v[pCO2].c[0] = v[pCO2].c[1];
            }
        }
        diagnostics_validate_tracer_state("Biogeo:exit", t);
        return;
    }
    


    /* ================================================================
     * MAIN LOOP START
     * ================================================================ */

    // Special Optimization: Stage 2 (SPM + P Adsorption) ONLY
    if (enable_biogeochemical_reactions == BIOGEO_P_ADSORPTION_ONLY) {
        const double dt = (double)DELTI;
        for (int i = 1; i <= M; i++) {
             compute_phosphorus_adsorption_step(i, dt);
        }
        diagnostics_validate_tracer_state("Biogeo:exit_P_only", t);
        return;
    }

    const double dt = (double)DELTI;
    const double I0_val = I0(t);

    // ---------------------------------------------------------------
    // Discharge-driven non-algal turbidity:
    //   KD_flow = k_turb_Q × max(0, Q_ups/Q_ref - 1)^alpha_turb
    // Q_turb_ref auto-computed from forcing mean if not set.
    // Ref: Uncles et al. (2002); Fettweis et al. (2006)
    static double Q_turb_ref_effective = 0.0;
    static int Q_turb_ref_initialized = 0;
    if (!Q_turb_ref_initialized) {
        if (Q_turb_ref > 0.0) {
            Q_turb_ref_effective = Q_turb_ref;
        } else if (forcingData[FORCING_DISCHARGE].dataSize > 0 &&
                   forcingData[FORCING_DISCHARGE].data) {
            double sum_Q = 0.0;
            int n = forcingData[FORCING_DISCHARGE].dataSize;
            for (int k = 0; k < n; k++) {
                sum_Q += fabs(forcingData[FORCING_DISCHARGE].data[k]);
            }
            Q_turb_ref_effective = sum_Q / n;
            if (Q_turb_ref_effective < 1.0) Q_turb_ref_effective = 1.0;
            if (debug_level >= DEBUG_LEVEL_ESSENTIAL) {
                printf("[INFO] Q_turb_ref auto-computed from forcing data: %.1f m3/s\n",
                       Q_turb_ref_effective);
            }
        } else {
            Q_turb_ref_effective = 1.0;  // Minimal fallback (effectively disables)
        }
        Q_turb_ref_initialized = 1;
    }

    double KD_flow = 0.0;
    if (k_turb_Q > 0.0 && Q_turb_ref_effective > 0.0) {
        double Q_ups = fabs(Discharge_ups(t));
        double q_excess = Q_ups / Q_turb_ref_effective - 1.0;
        if (q_excess > 0.0) {
            KD_flow = k_turb_Q * pow(q_excess, alpha_turb);
        }
    }
    
    /* ================================================================
     * PRE-COMPUTE UPSTREAM VOLUME AND NET FRESHWATER DISCHARGE
     * C-GEM grid: i=1 is Downstream (Mouth), i=M is Upstream (River).
     * Therefore, the volume "upstream" of cell i is the integral
     * from M down to i. The net discharge through cell i is the 
     * upstream boundary discharge plus any tributaries joining >= i.
     * ================================================================ */
    double vol_upstream_array[MAXM + 1] = {0.0};
    double Q_net_inflow_array[MAXM + 1] = {0.0};
    
    double accum_vol = 0.0;
    double current_Q_inflow = fabs(Discharge_ups(t)); 
    
    for (int i = M; i >= 1; i--) {
        if (freeArea[i] > 0.0) {
            accum_vol += freeArea[i] * (double)DELXI;
        }
        vol_upstream_array[i] = accum_vol;
        
        if (tributaryEnabled && numTributaries > 0 && tributaries) {
            for (int trib = 0; trib < numTributaries; ++trib) {
                if (tributaries[trib].cellIndex == i) {
                    double q = calculatePhysicsBasedTributaryDischarge(t, trib);
                    if (q > 0.0) current_Q_inflow += q;
                }
            }
        }
        Q_net_inflow_array[i] = current_Q_inflow;
    }

    /* ================================================================
     * MAIN SPATIAL LOOP - All cells from 1 to M
     * ================================================================ */
    for (int i = 1; i <= M; i++) {
        double sal = v[Sal].c[i];
        double depth = waterDepth[i];
        
        // Skip invalid cells
        if (depth <= 0.1) continue;
        
        // Ensure minimum salinity for stability
        double sal_safe = (sal < 1e-10) ? 0.0 : sal;
        
        /* ============================================================
         * 1. LIGHT ATTENUATION AND PRIMARY PRODUCTION
         * ============================================================ */
        // Light attenuation: 6-component process-based formulation
        //   (1) kbg: background mineral/water absorption [m⁻¹]
        //   (2) kspm × SPM: mass-specific particulate attenuation (Kirk 2010)
        //   (3) kCDOM × DOC: chromophoric dissolved organic matter (Fichot & Benner 2012)
        //       DOC proxy: TOC × 12.01/1000 converts mmolC/m³ → mgC/L
        //   (4) KD_Phy × Phy: phytoplankton self-shading (Cloern 1987)
        //   (5) KD_flow: discharge-driven turbidity (computed above)
        //   (6) KD_resus: depth-dependent resuspension turbidity
        //       (Dyer 1997; Winterwerp 2006): KD_resus = k_resus / depth
        //
        // Spatially varying background light attenuation
        // Two complementary mechanisms contribute to the longitudinal
        // turbidity gradient observed in tropical estuaries:
        //
        // (a) Distance-based terrigenous CDOM/mineral turbidity:
        //     Upstream sections carry terrigenous humic/fulvic CDOM and
        //     fine clay particles from laterite/ferralsol catchment soils.
        //     The transition from clear marine-influenced water to turbid
        //     freshwater is parameterized as a smooth logistic function
        //     of distance from the estuary mouth, consistent with observed
        //     KD profiles from satellite remote sensing (Doxaran et al. 2009).
        //     kbg_upstream : excess upstream attenuation [m⁻¹]
        //     kbg_transition_km : distance of turbidity transition [km]
        //     kbg_blend_km : half-width of transition zone [km]
        //
        // (b) Salinity-dependent CDOM dilution (optional, additive):
        //     Refs: Blough & Del Vecchio (2002), Fichot & Benner (2012)
        double kbg_local = kbg;
        // (a) Distance-based upstream turbidity
        if (kbg_upstream > 0.0 && kbg_transition_km > 0.0) {
            double dist_km = (double)(i - 1) * DELXI / 1000.0;
            double blend = (kbg_blend_km > 0.0) ? kbg_blend_km : 5.0;
            double f_upstream = 0.5 * (1.0 + tanh((dist_km - kbg_transition_km) / blend));
            kbg_local += kbg_upstream * f_upstream;
        }
        // (b) Salinity-dependent freshwater CDOM (additive)
        if (kbg_fresh > 0.0 && S_cdom_width > 0.0) {
            double sal_local = (sal_safe > 0.0) ? sal_safe : 0.0;
            double cdom_frac = 0.5 * (1.0 + tanh((S_cdom_threshold - sal_local) / S_cdom_width));
            kbg_local += kbg_fresh * cdom_frac;
        }
        // Depth-dependent resuspension turbidity
        double KD_resus = 0.0;
        if (k_resus > 0.0 && depth > 0.5) {
            KD_resus = k_resus / depth;
        }
        const double KD_turbidity = kbg_local + kspm * (1000.0 * v[SPM].c[i]) + KD_flow + KD_resus;
        double KD_cdom = kCDOM * (v[TOC].c[i] * 12.01 / 1000.0);  // TOC as DOC proxy
        if (!isfinite(KD_cdom) || KD_cdom < 0.0) KD_cdom = 0.0;
        double KD = KD_turbidity + KD_cdom + KD_Phy * (v[Phy1].c[i] + v[Phy2].c[i]);
        
        // Nutrient limitation for Phy1 (diatoms - includes Si)
        double N_total = v[NO3].c[i] + v[NH4].c[i];
        double nlim_1 = (v[Si].c[i] / (v[Si].c[i] + KSi[Phy1])) *
                        (N_total / (N_total + KN[Phy1])) *
                        (v[PO4].c[i] / (v[PO4].c[i] + KPO4[Phy1]));
        
        // Nutrient limitation for Phy2 (non-diatoms - no Si limitation)
        double nlim_2 = (N_total / (N_total + KN[Phy2])) *
                        (v[PO4].c[i] / (v[PO4].c[i] + KPO4[Phy2]));

        // Mechanistic CO2 limitation (disabled by default: enable_co2_limitation=0).
        // In turbid tropical estuaries with pCO2 >> K_CO2, dissolved CO2 rarely
        // limits phytoplankton growth.  Enable for systems with lower pCO2.
        double co2lim = 1.0;
        if (enable_co2_limitation && (K_CO2_Phy > 0.0) && isfinite(K_CO2_Phy)) {
            double co2 = v[CO2].c[i];
            if (!isfinite(co2) || co2 < 0.0) co2 = 0.0;
            co2lim = co2 / (co2 + K_CO2_Phy);
            if (!isfinite(co2lim) || co2lim < 0.0) co2lim = 0.0;
            if (co2lim > 1.0) co2lim = 1.0;
        }
        
        // N-switch for NH4 preference (tunable K_NH4_switch)
        double nswitch = v[NH4].c[i] / (K_NH4_switch + v[NH4].c[i]);
        
        // Temperature-dependent maximum photosynthetic rate (using existing functions)
        double Pbmax_1 = Pbmax(t, Phy1);
        double Pbmax_2 = Pbmax(t, Phy2);
        
        // Light limitation using clean Fortran formulation
        double lightlim_1 = light_limitation(alpha[Phy1], Pbmax_1, I0_val, KD, depth);
        double lightlim_2 = light_limitation(alpha[Phy2], Pbmax_2, I0_val, KD, depth);
        
        // Store diagnostic limitation factors (Phy1 / diatoms as primary diagnostic)
        diag_fN[i]  = N_total / (N_total + KN[Phy1]);
        diag_fP[i]  = v[PO4].c[i] / (v[PO4].c[i] + KPO4[Phy1]);
        diag_fSi[i] = v[Si].c[i] / (v[Si].c[i] + KSi[Phy1]);
        diag_fI[i]  = (depth > 0.0) ? lightlim_1 / depth : 0.0;  // depth-averaged [0,1]
        diag_KD_total[i] = KD;

        // Gross Primary Production (GPP) [mmol C m⁻² s⁻¹]
        // CO2 limitation term (co2lim) — currently always 1.0 (feature disabled).
        // Salinity effects handled via sal_death mortality (Arndt et al. 2007),
        // not growth inhibition — freshwater Phy are killed by high salinity,
        // not prevented from growing.
        double GPP_1 = Pbmax_1 * v[Phy1].c[i] * nlim_1 * lightlim_1 * co2lim;
        double GPP_2 = Pbmax_2 * v[Phy2].c[i] * nlim_2 * lightlim_2 * co2lim;
        
        // Temperature-dependent maintenance rate (using existing function)
        double kmaint_1 = kmaintenance(t, Phy1);
        double kmaint_2 = kmaintenance(t, Phy2);
        
        // Net Primary Production split by N source [mmol C m⁻³ s⁻¹]
        // NOTE: Maintenance is subtracted ONCE from total NPP, not from each N-source pathway
        // Original formulation: NPP = GPP * (1-excr) * (1-growth) / depth - maintenance
        double gross_growth_1 = (GPP_1 / depth) * (1.0 - kexcr[Phy1]) * (1.0 - kgrowth[Phy1]);
        double gross_growth_2 = (GPP_2 / depth) * (1.0 - kexcr[Phy2]) * (1.0 - kgrowth[Phy2]);
        
        // Split gross growth by N source (NO3 vs NH4)
        double NPP_NO3_1 = (1.0 - nswitch) * gross_growth_1;
        double NPP_NH4_1 = nswitch * gross_growth_1;
        double NPP_NO3_2 = (1.0 - nswitch) * gross_growth_2;
        double NPP_NH4_2 = nswitch * gross_growth_2;
        
        // Maintenance respiration (subtracted separately, not from each pathway)
        double maint_1 = kmaint_1 * v[Phy1].c[i];
        double maint_2 = kmaint_2 * v[Phy2].c[i];
        
        // Total NPP
        double NPP_NO3_total = NPP_NO3_1 + NPP_NO3_2;
        double NPP_NH4_total = NPP_NH4_1 + NPP_NH4_2;
        double NPP_total = NPP_NO3_total + NPP_NH4_total;
        
        /* ============================================================
         * PHYTOPLANKTON MORTALITY — Standard NPZD formulation
         * ============================================================
         * Loss terms:
         *   (1) Linear mortality: temperature-dependent natural death (Eppley 1972)
         *   (2) Quadratic closure: implicit zooplankton grazing + viral lysis
         *       (Fennel et al. 2006 JGR; Edwards & Yool 2000 JPR)
         *   (3) Salinity stress: osmotic mortality for freshwater taxa in
         *       brackish/saline zones (Garnier et al., RIVE; Billen et al. 2001)
         *       f(S) = S^2 / (S^2 + S_tol^2): 0 in freshwater, ~1 in saline
         * ============================================================ */
        const double kmort_base_1 = kmort(t, Phy1);
        const double kmort_base_2 = kmort(t, Phy2);

        // Linear mortality (natural death) → products stay in water column (TOC)
        double linear_death_1 = kmort_base_1 * v[Phy1].c[i];
        double linear_death_2 = kmort_base_2 * v[Phy2].c[i];

        // Salinity stress mortality (freshwater taxa die in saline water)
        double sal_death_1 = 0.0, sal_death_2 = 0.0;
        if (kmort_sal > 0.0 && sal_safe > 0.1) {
            double f_sal_1 = (sal_safe * sal_safe) / (sal_safe * sal_safe + sal_tol_Phy1 * sal_tol_Phy1);
            double f_sal_2 = (sal_safe * sal_safe) / (sal_safe * sal_safe + sal_tol_Phy2 * sal_tol_Phy2);
            sal_death_1 = kmort_sal * f_sal_1 * v[Phy1].c[i];
            sal_death_2 = kmort_sal * f_sal_2 * v[Phy2].c[i];
        }

        // Quadratic mortality closure → products exported as sinking fecal pellets
        // Remineralized in benthos, returned via SOD/benthic recycling terms.
        double quadratic_death_1 = kmort2_Phy1 * v[Phy1].c[i] * v[Phy1].c[i];
        double quadratic_death_2 = kmort2_Phy2 * v[Phy2].c[i] * v[Phy2].c[i];

        double death_1 = linear_death_1 + quadratic_death_1 + sal_death_1;
        double death_2 = linear_death_2 + quadratic_death_2 + sal_death_2;
        double death_total = death_1 + death_2;
        double linear_death_total = linear_death_1 + linear_death_2 + sal_death_1 + sal_death_2;
        
        // Si consumption (diatoms only). Computed after any rate limiting so that
        // Si uptake remains consistent with the final (possibly scaled) diatom NPP.
        double Si_cons = 0.0;
        
        /* ============================================================
         * 2. SECONDARY BIOGEOCHEMICAL REACTIONS
         * ============================================================ */
        
        // Aerobic degradation (using existing Fhetox function)
        double R_adegrad = 0.0;
        if (v[TOC].c[i] > 1e-10 && v[O2].c[i] > 1e-10) {
            R_adegrad = Fhetox(t) * 
                        (v[TOC].c[i] / (v[TOC].c[i] + KTOC)) *
                        (v[O2].c[i] / (v[O2].c[i] + KO2_ox));
        }
        
        // Denitrification (using existing Fhetden function)
        double R_denit = 0.0;
        if (v[TOC].c[i] > 1e-10 && v[NO3].c[i] > 1e-10) {
            R_denit = Fhetden(t) *
                      (v[TOC].c[i] / (v[TOC].c[i] + KTOC)) *
                      (KinO2 / (v[O2].c[i] + KinO2)) *
                      (v[NO3].c[i] / (v[NO3].c[i] + KNO3));
        }
        
        // Nitrification (using existing Fnit function)
        double R_nitrif = 0.0;
        if (v[NH4].c[i] > 1e-10 && v[O2].c[i] > 1e-10) {
            R_nitrif = Fnit(t) *
                       (v[O2].c[i] / (v[O2].c[i] + KO2_nit)) *
                       (v[NH4].c[i] / (v[NH4].c[i] + KNH4));
        }

        // Fast labile-C oxygen demand: dual-rate TOC degradation (Arndt et al. 2011)
        // Represents the labile fraction of bulk TOC that degrades faster than the
        // refractory pool (kox). Applied globally with Monod O2/TOC kinetics.
        double R_cbod_fast = 0.0;
        if (kcbod_fast > 0.0 && v[TOC].c[i] > 1e-10 && v[O2].c[i] > 1e-10) {
            const double f_toc_fast = v[TOC].c[i] / (v[TOC].c[i] + KTOC_fast);
            const double f_o2_fast = v[O2].c[i] / (v[O2].c[i] + KO2_cbod);
            R_cbod_fast = kcbod_fast * f_toc_fast * f_o2_fast * v[TOC].c[i];
        }

        double R_om_ox = R_adegrad + R_cbod_fast;

        // (R_ww_bod removed: O2 demand must come from state-variable reactions for mass conservation)
        
        // O2 air-water exchange (using piston velocity)
        double PisVel = piston_velocity(i, t);
        double R_O2_ex = (PisVel / depth) * (O2sat(t, i) - v[O2].c[i]);
        
        /* ============================================================
         * 3. pH AND CARBONATE CHEMISTRY (Robust DIC/TA -> pH solver)
         * ============================================================ */

        double pH_now = v[PH].c[i];
        double CO2aq_eq = v[CO2].c[i];
        double pCO2_now = v[pCO2].c[i];
        int carbonate_ok = 0;
        int carbonate_ok_for_flux = 0;

        // Performance switch: Stage 3 calibration typically does not target carbonate.
        // In normal (non-calibration) runs, carbonate chemistry is part of BIOGEO_FULL
        // regardless of the diagnostics flag (tests rely on CO2 exchange behavior).
        // During calibration, we allow disabling carbonate to speed optimization.
        const int do_carbonate = (calibration_mode > 0) ? enable_carbonate_diagnostics : 1;

        if (do_carbonate) {
            carbonate_ok = carbonate_solve_dic_ta(v[DIC].c[i], v[AT].c[i], v[NH4].c[i], v[NO3].c[i], sal, t, i, &pH_now, &CO2aq_eq, &pCO2_now);
            carbonate_ok_for_flux = carbonate_ok;
            if (!carbonate_ok) {
                carbonate_solver_failures_day++;
                if (carbonate_solver_failure_policy == CARBONATE_FAIL_FATAL) {
                    fprintf(stderr,
                            "FATAL: Carbonate solver failed at cell=%d t=%d (DIC=%.17g AT=%.17g Sal=%.17g)\n",
                            i, t, v[DIC].c[i], v[AT].c[i], sal);
                    exit(EXIT_FAILURE);
                }
                if (carbonate_solver_failure_policy == CARBONATE_FAIL_SKIP_CO2_EXCHANGE) {
                    carbonate_ok_for_flux = 0;
                }
                // Fall back to last-known diagnostics (do not abort; keep model running)
                pH_now = v[PH].c[i];
                CO2aq_eq = v[CO2].c[i];
                pCO2_now = v[pCO2].c[i];
            } else {
                v[PH].c[i] = pH_now;
                v[CO2].c[i] = CO2aq_eq;
                v[pCO2].c[i] = pCO2_now;
            }
        }

        // CO2 air-water exchange (0.913 factor for CO2 vs O2 Schmidt number ratio)
        double R_CO2_ex = 0.0;
        if (do_carbonate) {
            double Henry = KH(t, i);
            // Atmospheric pCO2 uses µatm (see assignBiogeochemicalRateConstants()); convert to atm via /1e6.
            const double CO2aq_eq_atm = Henry * (pCO2atmo / 1e6);
            if (carbonate_ok_for_flux || carbonate_solver_failure_policy == CARBONATE_FAIL_HOLD_LAST) {
                R_CO2_ex = -(PisVel * 0.913 / depth) * (CO2aq_eq - CO2aq_eq_atm);
            } else {
                // SKIP_CO2_EXCHANGE and carbonate solver failed: do not apply gas exchange this step.
                R_CO2_ex = 0.0;
            }
        }
        
        // DIC/TA reaction rates are computed after limiting so they remain
        // consistent with the final (applied) process rates.
        double R_DIC = 0.0;
        double R_TA = 0.0;

        /* ============================================================
         * 4b. CONSOLIDATED METABOLIC SAFEGUARDS (Budget-Preserving)
         * ============================================================
         * To prevent negative concentrations under aggressive calibration, we cap 
         * reaction rates based on available supply. Guards are ordered by 
         * dependency to ensure supply-chain integrity.
         */
        if (dt > 0.0) {
            // --- 1. TOC & O2 (Limiting Aerobic Degradation & Denitrification) ---
            // These consume TOC/O2 but provide supply for NH4/NO3/PO4.
            const double total_toc_demand = R_om_ox + R_denit;
            const double max_toc_consumption = (v[TOC].c[i] / dt) + death_total;
            if (total_toc_demand > 0.0 && isfinite(max_toc_consumption)) {
                if (max_toc_consumption <= 0.0) {
                    R_adegrad = 0.0; R_cbod_fast = 0.0; R_denit = 0.0;
                } else if (total_toc_demand > max_toc_consumption) {
                    const double scale = max_toc_consumption / total_toc_demand;
                    R_adegrad *= scale; R_cbod_fast *= scale; R_denit *= scale;
                }
            }
            R_om_ox = R_adegrad + R_cbod_fast;

            const double oxygen_supply_rate = (v[O2].c[i] / dt) + NPP_NH4_total + (138.0/106.0)*NPP_NO3_total + R_O2_ex;
            const double oxygen_demand_rate = R_om_ox + 2.0 * R_nitrif;
            if (oxygen_demand_rate > 0.0 && isfinite(oxygen_supply_rate)) {
                if (oxygen_supply_rate <= 0.0) {
                    R_adegrad = 0.0; R_cbod_fast = 0.0; R_nitrif = 0.0;
                } else if (oxygen_demand_rate > oxygen_supply_rate) {
                    const double scale = oxygen_supply_rate / oxygen_demand_rate;
                    R_adegrad *= scale; R_cbod_fast *= scale; R_nitrif *= scale;
                }
            }
            R_om_ox = R_adegrad + R_cbod_fast;

            // --- 2. AMMONIUM (Limiting Nitrification & NH4 uptake) ---
            // Supply = (NH4/dt) + mineralization from R_adegrad using CN_toc stoichiometry.
            // R_cbod_fast is decoupled from N production (C-rich labile fraction;
            // QUAL2E-style separation of BOD kinetics from organic-N kinetics).
            const double redn_toc_cap = 1.0 / CN_toc;  // N:C for bulk TOC
            const double nh4_supply_rate = (v[NH4].c[i] / dt) + redn_toc_cap * R_adegrad;
            const double nh4_demand_rate = redn * NPP_NH4_total + R_nitrif;
            if (nh4_demand_rate > 0.0 && isfinite(nh4_supply_rate)) {
                if (nh4_supply_rate <= 0.0) {
                    NPP_NH4_1 = 0.0; NPP_NH4_2 = 0.0; R_nitrif = 0.0;
                } else if (nh4_demand_rate > nh4_supply_rate) {
                    const double scale = nh4_supply_rate / nh4_demand_rate;
                    NPP_NH4_1 *= scale; NPP_NH4_2 *= scale; R_nitrif *= scale;
                }
                NPP_NH4_total = NPP_NH4_1 + NPP_NH4_2; // Update for final NPP total
            }

            // --- 3. NITRATE (Limiting Denitrification & NO3 uptake) ---
            const double denit_no3_coeff = (94.4 / 106.0);
            const double other_rate_no3 = (R_nitrif - redn * NPP_NO3_total); // Nitrif now capped
            const double max_denit_by_no3 = (denit_no3_coeff > 0.0) ? ((v[NO3].c[i] / dt) + other_rate_no3) / denit_no3_coeff : 0.0;
            if (R_denit > 0.0 && isfinite(max_denit_by_no3)) {
                if (max_denit_by_no3 <= 0.0) {
                    R_denit = 0.0;
                } else if (R_denit > max_denit_by_no3) {
                    R_denit = max_denit_by_no3;
                }
            }
            // Also cap NO3 uptake to avoid independent NO3 crash
            const double no3_supply_for_uptake = (v[NO3].c[i] / (redn * dt)) + (R_nitrif / redn) - (denit_no3_coeff * R_denit / redn);
            if (NPP_NO3_total > 0.0 && isfinite(no3_supply_for_uptake)) {
                if (no3_supply_for_uptake <= 0.0) {
                    NPP_NO3_1 = 0.0; NPP_NO3_2 = 0.0;
                } else if (NPP_NO3_total > no3_supply_for_uptake) {
                    const double scale = no3_supply_for_uptake / NPP_NO3_total;
                    NPP_NO3_1 *= scale; NPP_NO3_2 *= scale;
                }
                NPP_NO3_total = NPP_NO3_1 + NPP_NO3_2;
            }

            // Sync total NPP after all Nitrogen/Carbon adjustments
            NPP_total = NPP_NO3_total + NPP_NH4_total;

            // --- 4. PHOSPHORUS (Final NPP scaling) ---
            const double po4_supply_rate = (v[PO4].c[i] / (redp * dt)) + R_adegrad + R_denit;
            if (NPP_total > 0.0 && isfinite(po4_supply_rate)) {
                if (po4_supply_rate <= 0.0) {
                    NPP_NO3_1 = 0.0; NPP_NO3_2 = 0.0; NPP_NH4_1 = 0.0; NPP_NH4_2 = 0.0;
                } else if (NPP_total > po4_supply_rate) {
                    const double scale = po4_supply_rate / NPP_total;
                    NPP_NO3_1 *= scale; NPP_NO3_2 *= scale; NPP_NH4_1 *= scale; NPP_NH4_2 *= scale;
                }
                // Final totals re-sync
                NPP_NO3_total = NPP_NO3_1 + NPP_NO3_2;
                NPP_NH4_total = NPP_NH4_1 + NPP_NH4_2;
                NPP_total = NPP_NO3_total + NPP_NH4_total;
            }

            // --- 5. SILICA (Diatom specific scaling) ---
            const double cur_NPP_Phy1 = NPP_NH4_1 + NPP_NO3_1;
            const double max_NPP_Phy1_by_si = v[Si].c[i] / (redsi * dt);
            if (cur_NPP_Phy1 > 0.0 && isfinite(max_NPP_Phy1_by_si)) {
                if (max_NPP_Phy1_by_si <= 0.0) {
                    NPP_NH4_1 = 0.0; NPP_NO3_1 = 0.0;
                } else if (cur_NPP_Phy1 > max_NPP_Phy1_by_si) {
                    const double scale = max_NPP_Phy1_by_si / cur_NPP_Phy1;
                    NPP_NH4_1 *= scale; NPP_NO3_1 *= scale;
                }
                // Last updates to totals
                NPP_NO3_total = NPP_NO3_1 + NPP_NO3_2;
                NPP_NH4_total = NPP_NH4_1 + NPP_NH4_2;
                NPP_total = NPP_NO3_total + NPP_NH4_total;
            }
        }

        // Update derived Si uptake based on final diatom NPP.
        Si_cons = -redsi * (NPP_NH4_1 + NPP_NO3_1);

        /* ============================================================
         * 4. SEDIMENT DYNAMICS
         * ============================================================ */
        compute_sediment_fluxes(i);
        
        /* ============================================================
         * 5. UPDATE CONCENTRATIONS [mmol m⁻³]
         * ============================================================ */
        
        // Phytoplankton Settling (RIVE alignment)
        // Loss from water column: - (ws_phy / H) * C
        // Prevents unrealistic biomass accumulation in upstream sections
        double H_safe = (waterDepth[i] > 0.5) ? waterDepth[i] : 1.0;
        double settle_rate = ws_phy / H_safe;

        // Benthic grazing — depth-dependent phytoplankton removal
        // Cloern (1982, Ecol. Monogr.): benthic filter-feeders access the full
        // water column in shallow sections but only the bottom layer in deep
        // channels. Loss rate = k_benth_graze / H creates differential removal:
        //   shallow upstream (H=10m): strong benthic control
        //   deep urban (H=22m): weak benthic control → bloom persists
        // Lucas et al. (1999, Estuaries): 0.5–10 m/day clearance in estuaries.
        double benth_graze_rate = k_benth_graze / H_safe;

        // Phytoplankton: budget-preserving update using Patankar scheme for sinks
        // (Ensures positivity even if rates are high relative to dt)
        // Patankar sink includes quadratic mortality linearized at current concentration
        // and salinity stress mortality rate (0 in freshwater, increasing in saline water)
        double sal_rate_1 = 0.0, sal_rate_2 = 0.0;
        if (kmort_sal > 0.0 && sal_safe > 0.1) {
            sal_rate_1 = kmort_sal * (sal_safe * sal_safe) / (sal_safe * sal_safe + sal_tol_Phy1 * sal_tol_Phy1);
            sal_rate_2 = kmort_sal * (sal_safe * sal_safe) / (sal_safe * sal_safe + sal_tol_Phy2 * sal_tol_Phy2);
        }
        double total_sink_rate_1 = kmaint_1 + kmort_base_1 + settle_rate + benth_graze_rate
                                 + kmort2_Phy1 * v[Phy1].c[i] + sal_rate_1;
        double total_sink_rate_2 = kmaint_2 + kmort_base_2 + settle_rate + benth_graze_rate
                                 + kmort2_Phy2 * v[Phy2].c[i] + sal_rate_2;
        
        // C_new = (C_old + NPP * dt) / (1 + SinkRate * dt)
        v[Phy1].c[i] = (v[Phy1].c[i] + (NPP_NO3_1 + NPP_NH4_1) * dt) / (1.0 + total_sink_rate_1 * dt);
        v[Phy2].c[i] = (v[Phy2].c[i] + (NPP_NO3_2 + NPP_NH4_2) * dt) / (1.0 + total_sink_rate_2 * dt);
        
        // Silicate: diatom uptake + regeneration from dead diatom dissolution
        // Si_cons = uptake by NPP (negative), Si_regen = release from dead diatoms (positive)
        double Si_regen = redsi * linear_death_1 * ksi_diss;  // Biogenic Si dissolution (water-column only)
        v[Si].c[i] += (Si_cons + Si_regen) * dt;
        
        // Variable C:N stoichiometry for bulk TOC mineralization
        const double redn_toc_pho = 1.0 / CN_toc;
        const double redp_toc = (1.0 / 106.0) * (6.625 / CN_toc);
        
        // Phosphate: TOC mineralization uses CN_toc-scaled stoichiometry,
        // phytoplankton uptake uses Redfield C:P.
        v[PO4].c[i] += (redp_toc * (R_adegrad + R_denit) - redp * NPP_total) * dt;
        
        /* ============================================================
         * 5a. PHOSPHORUS ADSORPTION/DESORPTION (Langmuir Isotherm)
         * ============================================================
         */
        compute_phosphorus_adsorption_step(i, dt);
        
        // Sediment O2 demand (SOD): benthic respiration consumes O2
        // SOD_flux = SOD_rate [mmol/m²/s] * f(O2) * f(TOC) / depth [m] = [mmol/m³/s]
        // Mechanistic localization (no hard zoning):
        // - f(O2): Monod oxygen control (diagenetic demand limited under hypoxia)
        // - f(TOC): organic substrate availability proxy using water-column TOC
        //   (Soetaert-style benthic/OM coupling; keeps higher SOD in polluted reaches)
        double R_SOD = 0.0;
        if (SOD_rate > 0.0 && depth > 0.0 && v[O2].c[i] > 0.0) {
            const double f_o2 = v[O2].c[i] / (v[O2].c[i] + KO2_SOD);
            double f_toc = 1.0;
            if (KTOC > 0.0 && isfinite(KTOC) && isfinite(v[TOC].c[i]) && v[TOC].c[i] > 0.0) {
                f_toc = v[TOC].c[i] / (v[TOC].c[i] + KTOC);
            }
            if (!isfinite(f_toc) || f_toc < 0.0) f_toc = 0.0;
            if (f_toc > 1.0) f_toc = 1.0;
            R_SOD = SOD_rate * f_o2 * f_toc / depth;
        }

        // Final oxygen-availability cap (positivity-preserving):
        // Include ALL aerobic O2 sinks that act in this step, including SOD and
        // maintenance respiration. This avoids numerical O2 undershoot during
        // global calibration while remaining mechanistic (oxygen-limited aerobic rates).
        if (dt > 0.0) {
            const double o2_sources = NPP_NH4_total + (138.0 / 106.0) * NPP_NO3_total + R_O2_ex;
            const double o2_sinks = R_om_ox + 2.0 * R_nitrif + R_SOD + (maint_1 + maint_2);
            const double o2_available = o2_sources + (v[O2].c[i] / dt);

            if (o2_sinks > 0.0 && isfinite(o2_available)) {
                if (o2_available <= 0.0) {
                    R_adegrad = 0.0;
                    R_cbod_fast = 0.0;
                    R_nitrif = 0.0;
                    R_SOD = 0.0;
                    maint_1 = 0.0;
                    maint_2 = 0.0;
                } else if (o2_sinks > o2_available) {
                    const double scale = o2_available / o2_sinks;
                    R_adegrad *= scale;
                    R_cbod_fast *= scale;
                    R_nitrif *= scale;
                    R_SOD *= scale;
                    maint_1 *= scale;
                    maint_2 *= scale;
                }
                R_om_ox = R_adegrad + R_cbod_fast;
            }
        }

        // Benthic remineralization coupled to SOD oxygen demand.
        // Uses CN_toc stoichiometry (redn_toc_pho / redp_toc defined above).
        double R_benthic_DIC = R_SOD;
        double R_benthic_NH4 = redn_toc_pho * R_SOD;
        double R_benthic_PO4 = redp_toc * R_SOD;
        
        // Oxygen: respiration + nitrification + photosynthesis + gas exchange - SOD - maintenance - wastewater BOD
        v[O2].c[i] += (-R_om_ox + NPP_NH4_total + 
                      (138.0 / 106.0) * NPP_NO3_total - 
                      2.0 * R_nitrif + R_O2_ex - R_SOD - (maint_1 + maint_2)) * dt;
        
        // TOC: linear mortality - degradation - denitrification
        // Quadratic mortality products are exported (sinking fecal pellets, Fennel 2006)
        // and remineralized in the sediment (returned via SOD/benthic recycling)
        v[TOC].c[i] += (-R_om_ox - R_denit + linear_death_total) * dt;
        
        // Ammonium: mineralization uses CN_toc (bulk TOC stoichiometry),
        // while phytoplankton uptake uses Redfield CN (Phy stoichiometry).
        // R_cbod_fast oxidizes C-rich labile organics (carbohydrates) without
        // NH4 release — consistent with QUAL2E BOD/org-N separation.
        v[NH4].c[i] += (redn_toc_pho * R_adegrad - redn * NPP_NH4_total - R_nitrif + R_benthic_NH4) * dt;
        
        // Phosphate: benthic recycling (RIVE-style, Garnier et al.)
        // Sediment mineralization releases P proportional to SOD via Redfield C:P ratio.
        v[PO4].c[i] += R_benthic_PO4 * dt;
        
        // Nitrate: nitrification - denitrification - uptake
        v[NO3].c[i] += (-(94.4 / 106.0) * R_denit + R_nitrif - 
                       redn * NPP_NO3_total) * dt;

         // DIC reaction rate (CO2 exchange + degradation + denitrification + maintenance - photosynthesis)
         R_DIC = R_CO2_ex + R_om_ox + R_denit + (maint_1 + maint_2) - NPP_total;

         // TA reaction rate (Soetaert et al. 2007, Table 2)
         // NH4 production from TOC mineralization uses CN_toc stoichiometry;
         // phytoplankton uptake uses Redfield.
         R_TA = redn_toc_pho * (R_adegrad) +                              // NH4 from TOC mineralization
             (15.0 / 106.0) * (maint_1 + maint_2) +                   // NH4 from maintenance (Phy → Redfield)
             (93.4 / 106.0) * R_denit -                                // NO3 consumption
             2.0 * R_nitrif +                                          // H+ production
             (-15.0 / 106.0) * NPP_NH4_total +                        // NH4 uptake
             (17.0 / 106.0) * NPP_NO3_total;                          // NO3 uptake
        
        // DIC and Alkalinity
        v[DIC].c[i] += (R_DIC + R_benthic_DIC) * dt;
        v[AT].c[i] += R_TA * dt;

        // Recompute carbonate diagnostics *after* updating DIC/AT (diagnostic treatment).
        // This avoids drift/inconsistency from treating CO2/pCO2 as prognostic variables.
        if (do_carbonate) {
            double pH2 = v[PH].c[i];
            double CO2_2 = v[CO2].c[i];
            double pCO2_2 = v[pCO2].c[i];
            int ok2 = carbonate_solve_dic_ta(v[DIC].c[i], v[AT].c[i], v[NH4].c[i], v[NO3].c[i], sal_safe, t, i, &pH2, &CO2_2, &pCO2_2);
            if (!ok2) {
                carbonate_solver_failures_day++;
            } else {
                v[PH].c[i] = pH2;
                v[CO2].c[i] = CO2_2;
                v[pCO2].c[i] = pCO2_2;
            }

            if (!isfinite(v[pCO2].c[i]) || v[pCO2].c[i] < 0.0 ||
                !isfinite(v[PH].c[i]) || v[PH].c[i] < 0.0 ||
                !isfinite(v[CO2].c[i]) || v[CO2].c[i] < 0.0) {
                carbonate_nonfinite_day++;
                diagnostics_fatal_tracer_error("Biogeo:carbonate_diag", t, pCO2, i, v[pCO2].c[i], "nonfinite_or_negative");
                return;
            }

            // Broad diagnostic flag (do not clamp; only count extremes for review)
            if (v[pCO2].c[i] > 1e6 || v[PH].c[i] < 4.0 || v[PH].c[i] > 10.5) {
                carbonate_extreme_day++;
            }
        }
        
        // NOTE: SPM erosion/deposition is handled by updateSuspendedSediment() in main loop
        // Do not duplicate SPM calculation here to avoid double-counting
        
        /* ============================================================
         * 6. NON-NEGATIVITY CONSTRAINTS
         * ============================================================ */
        // Enforced globally by diagnostics_validate_tracer_state() (hard-abort).
        
        /* ============================================================
         * 7. STORE REACTION RATES FOR OUTPUT
         * ============================================================ */
        adegrad[i] = R_adegrad;
        denit[i] = R_denit;
        nitrif[i] = R_nitrif;
        o2air[i] = R_O2_ex;
        co2air[i] = R_CO2_ex;
        reactionDIC[i] = R_DIC;
        reactionTA[i] = R_TA;
        NPP[i] = NPP_total;
        NPP_NO3_tot[i] = NPP_NO3_total;
        NPP_NH4_tot[i] = NPP_NH4_total;
        phydeath_tot[i] = death_total;
        
    }  // End spatial loop
    
    /* ================================================================
     * OUTPUT REACTION RATES (only at output intervals)
     * ================================================================ */
    if (t % (TS * DELTI) == 0 && t >= WARMUP) {
        Rates(NPP_NO3_tot, "Reaction_NPP_NO3", t);
        Rates(NPP_NH4_tot, "Reaction_NPP_NH4", t);
        Rates(phydeath_tot, "Reaction_phydeath", t);
        Rates(adegrad, "Reaction_adegradation", t);
        Rates(denit, "Reaction_denitrification", t);
        Rates(nitrif, "Reaction_nitrification", t);
        Rates(o2air, "Reaction_O2_exchange", t);
        Rates(co2air, "Reaction_CO2_exchange", t);
        Rates(reactionDIC, "Reaction_DIC", t);
        Rates(reactionTA, "Reaction_TA", t);
        // Diagnostic limitation factors (dimensionless 0-1) and KD (m⁻¹)
        Rates(diag_fN, "Diag_fN", t);
        Rates(diag_fP, "Diag_fP", t);
        Rates(diag_fSi, "Diag_fSi", t);
        Rates(diag_fI, "Diag_fI", t);
        Rates(diag_KD_total, "Diag_KD", t);
    }

    // =====================================================================
    // STAGGERED GRID CONSISTENCY: Post-Biogeo Smoothing
    // =====================================================================
    for (int s = 0; s < CHEM_COUNT; s++) {
        if (v[s].env != 1) continue;
        if (s == pCO2 || s == PH || s == CO2) continue; 
        
        for (int i = 2; i <= M2; i += 2) {
            v[s].c[i] = 0.5 * (v[s].c[i - 1] + v[s].c[i + 1]);
        }
        if (M % 2 == 0) v[s].c[M] = v[s].c[M - 1];
    }

    diagnostics_validate_tracer_state("Biogeo:exit", t);
}
