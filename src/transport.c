/**
 * @file transport.c
 * @brief Transport of chemical species (Advection-Dispersion)
 *        Method: Geometry-Driven Savenije (2012) with Recursive Dispersion
 * @author An Nguyen
 * @date Last modified: 12/2025
 * 
 * @references
 *   - Savenije, H.H.G. (2012). Salinity and Tides in Alluvial Estuaries. 2nd Ed.
 *   - Gisen, J.I.A., Savenije, H.H.G., Nijzink, R.C. (2015). Revised predictive equations
 *     for salt intrusion modelling in estuaries. HESS 19, 2791-2803.
 *   - Fischer, H.B., List, E.J., Koh, R.C.Y., et al. (1979). Mixing in Inland and
 *     Coastal Waters. Academic Press, New York.
 *   - Rutherford, J.C. (1994). River Mixing. Wiley, Chichester.
 */

#include "define.h"
#include "variables.h"
#include "diagnostics.h"
#include "file.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <ctype.h>
#include <errno.h>

// --- External Dependencies ---
extern double Discharge_ups(int t);
extern double Discharge(int t, int i, double *trib_discharge);
extern double calculatePhysicsBasedTributaryDischarge(int t, int tributary_index);
extern double interpolateInputData(int t, double* timeArray, double* dataArray, int dataSize, bool isHourly, const char* dataName, double defaultValue);

// --- Global Arrays from Hydrodynamics ---
extern double velocity[MAXM + 1];
extern double disp[MAXM + 1];
extern double width[MAXM + 1];
extern double velocity_min[MAXM + 1];
extern double velocity_max[MAXM + 1];
extern double LC1;      // Convergence length segment 1 [m]
extern double LC2;      // Convergence length segment 2 [m]
extern int index_2;     // Segment transition cell index
extern double C_VDB;    // Van der Burgh K scaling factor (calibration parameter)
extern double D0_CORRECTION; // Mouth dispersion correction factor (calibration parameter)


// Forward declaration
static double compute_velocity_amplitude(int cell_idx);
static int tributary_continuity_index(int cell);
static int upstream_tracer_index(void);



// ============================================================================
// PHYSICAL CONSTANTS (Literature-based, fixed values)
// ============================================================================

// Dispersion bounds are now SCALE-DEPENDENT (Gisen 2015 / Savenije 2012)
// D_min = 10% of mouth dispersion D0 (prevents unrealistic low values)
// D_max = 5× mouth dispersion D0 (prevents unrealistic high values)
// These are computed dynamically in Dispcoef() based on D0
static double DISPERSION_MIN_DYNAMIC = 50.0;   // Will be updated to 0.1 × D0
static double DISPERSION_MAX_DYNAMIC = 2000.0; // Will be updated to 5.0 × D0

// Absolute safety bounds (prevents numerical issues regardless of scale)
// Now configurable via params.txt (MAINT-001A fix)

// TVD flux limiter threshold (numerical, not physical)
static const double TVD_FLUX_THRESHOLD = 1.0e-35;  // Prevent division by zero in flux limiter

// Runtime diagnostics for reviewer package export
static int transport_time_context = 0;
static long tvd_calls_total = 0;
static long tvd_calls_post_warmup = 0;
static long tvd_substep_activation_events = 0;
static long tvd_substep_activation_events_post_warmup = 0;
static int tvd_max_substeps_observed = 1;
static long limiter_face_activation_events = 0;
static long limiter_face_activation_events_post_warmup = 0;
static long diffusion_substep_activation_events = 0;
static long diffusion_substep_activation_events_post_warmup = 0;

// --- Mouth Geometry Reference (initialized once) ---
static int mouth_geometry_initialized = 0;
static double mouth_width_ref = 1.0;
static double mouth_depth_ref = 5.0;
static double mouth_area_ref = 5.0;

/**
 * @brief Initialize mouth geometry reference values (once per simulation)
 */
static void ensure_mouth_geometry(void)
{
    if (mouth_geometry_initialized) return;

    // Width: use cell 1, fallback 100m
    mouth_width_ref = (width[1] > 1.0 && isfinite(width[1])) ? width[1] : 100.0;

    // Depth: prefer riverbed_depth, then waterDepth, fallback 5m
    double h = riverbed_depth[1];
    if (!(h > 0.5) || !isfinite(h)) h = waterDepth[1];
    if (!(h > 0.5) || !isfinite(h)) h = 5.0;
    mouth_depth_ref = h;

    // Area
    mouth_area_ref = mouth_width_ref * mouth_depth_ref;
    mouth_geometry_initialized = 1;
}

/**
 * @brief Sum all freshwater inputs (upstream + tributaries)
 * @note This is only used for logging and Canter-Cremers number.
 *       For cell-specific discharge, use Discharge(t, i, NULL) from bcforcing.c
 */
static double aggregate_freshwater_discharge(int t)
{
    double Q = fabs(Discharge_ups(t));
    if (tributaryEnabled && tributaries) {
        for (int i = 0; i < numTributaries; ++i) {
            double q = calculatePhysicsBasedTributaryDischarge(t, i);
            if (q > 0.0 && isfinite(q)) Q += q;
        }
    }
    return fmax(Q, 1.0);
}

/**
 * @brief Freshwater discharge magnitude used by Dispcoef recursion at section i.
 * @details
 * Uses mode-independent freshwater magnitude semantics so transport recursion is
 * physically consistent across tributary injection modes:
 *   - Main river contribution uses |Discharge_ups|.
 *   - Tributaries (positive inflow) are accumulated for confluences at or
 *     seaward of section i; with this model indexing (cell index increases
 *     upstream), that is tributary cellIndex >= i (same convention as Discharge()).
 */
static double dispersion_freshwater_magnitude(int t, int i)
{
    int section = i;
    if (section < 1) section = 1;
    if (section > M) section = M;

    double Q_mag = fabs(Discharge_ups(t));

    if (tributaryEnabled && tributaries && numTributaries > 0) {
        for (int j = 0; j < numTributaries; ++j) {
            if (tributaries[j].cellIndex >= section) {
                const double q = calculatePhysicsBasedTributaryDischarge(t, j);
                if (q > 0.0 && isfinite(q)) {
                    Q_mag += q;
                }
            }
        }
    }

    return fmax(Q_mag, 1e-9);
}

/**
 * @brief Canter-Cremers estuarine number (Savenije 2012 Eq. 5.6)
 * @details N = π × Q_freshwater / (h₀ × B₀)
 * 
 * The Canter-Cremers number is a dimensionless ratio comparing freshwater discharge
 * to tidal prism, characterizing estuary mixing behavior.
 * 
 * Note: This is the PURE literature formulation without any calibration factors.
 * Calibration is applied in the dispersion calculation via D0_CORRECTION.
 */
static double compute_canter_cremers_number(int t)
{
    ensure_mouth_geometry();
    double Q = aggregate_freshwater_discharge(t);
    // Savenije (2012) Eq. 5.6 - no calibration factor here, pure physics
    double N = (M_PI * Q) / fmax(mouth_depth_ref * mouth_width_ref, 1.0);
    return fmax(N, 1e-6);
}

/**
 * @brief Mouth dispersion coefficient D₀ (Savenije 2012 Eq. 5.72)
 * @details D₀ = 26 × h₀^1.5 × √(N × g) × D0_CORRECTION
 * 
 * D0_CORRECTION is a calibration parameter that accounts for site-specific
 * deviations from the empirical formula. Typical range: 0.5-2.0.
 */
static double compute_mouth_dispersion(int t)
{
    ensure_mouth_geometry();
    double N = compute_canter_cremers_number(t);
    // Savenije (2012) Eq. 5.72 with calibration via D0_CORRECTION
    double D0 = 26.0 * pow(mouth_depth_ref, 1.5) * sqrt(N * G) * D0_CORRECTION;
    
    // Discharge-dependent D0 damping (stratification correction):
    // In wet season (high N), estuary stratification reduces effective mixing,
    // attenuating D0 relative to the well-mixed Savenije (2012) prediction.
    // D0_eff = D0 × (N/N_ref)^(-K_DISCHARGE_SENSITIVITY/2)
    // This partially counteracts the √N scaling in the Savenije formula.
    if (K_DISCHARGE_SENSITIVITY > 0.0) {
        static double D0_N_ref = -1.0;
        static long   D0_step = 0;
        if (D0_N_ref < 0.0) {
            D0_N_ref = N;
        } else {
            double ema_alpha = 1.0 / 87600.0;
            D0_N_ref = D0_N_ref * (1.0 - ema_alpha) + N * ema_alpha;
        }
        D0_step++;
        if (D0_step > 24000 && D0_N_ref > 0.0) {
            // Half the K_DISCHARGE_SENSITIVITY for D0 damping (weaker than K correction)
            double damping = pow(N / D0_N_ref, -0.5 * K_DISCHARGE_SENSITIVITY);
            damping = fmin(fmax(damping, 0.3), 3.0);
            D0 *= damping;
        }
    }
    
    // Apply absolute safety bounds for D0 (dynamic bounds set later based on this D0)
    return fmin(fmax(D0, dispersion_abs_min), dispersion_abs_max);
}

/**
 * @brief Van der Burgh coefficient K (Gisen 2015 Eq. 40, calibrated)
 * @details K_predicted = 4.38 × h₀^0.36 × B₀^(-0.21) × Lc^(-0.14)
 *          K_calibrated = K_predicted × C_VDB (calibration scaling factor)
 * 
 * C_VDB is a calibration multiplier (NOT the Van der Burgh K itself).
 * The parameter controls dispersion decay rate along the estuary.
 * - C_VDB < 1.0: slower decay → more salt intrusion
 * - C_VDB > 1.0: faster decay → less salt intrusion
 * 
 * Literature K range: 0.2-0.8 (Savenije 2012 §4.4)
 * After C_VDB scaling, result is clamped to [0, 1]
 */
static double compute_van_der_burgh(double LC_local)
{
    ensure_mouth_geometry();
    // Gisen (2015) Eq. 40 - predictive formula
    double K_predicted = 4.38 * pow(mouth_depth_ref, 0.36) 
                              * pow(mouth_width_ref, -0.21) 
                              * pow(fmax(LC_local, 1.0), -0.14);
    // Apply calibration scaling factor C_VDB
    double K = K_predicted * C_VDB;
    // Clamp to physical bounds [0, 1] per Savenije (2012)
    return fmin(fmax(K, 0.0), 1.0);
}

/**
 * @brief Local cross-section convergence length: LCD = -Δx / ln(A[i]/A[i-1])
 */
static double compute_local_convergence_length(int idx, double fallback)
{
    if (idx <= 1 || idx > M) return fallback;
    
    double A_curr = totalArea[idx], A_prev = totalArea[idx - 1];
    if (!(A_curr > 0.0) || !(A_prev > 0.0)) return fallback;
    
    double ratio = A_curr / A_prev;
    if (fabs(ratio - 1.0) < 1e-6) return fallback;
    
    double lcd = -((double)DELXI) / log(ratio);
    return (lcd > 0.0 && isfinite(lcd)) ? lcd : fallback;
}

/**
 * @brief Tidal velocity amplitude (Savenije 2012 §4.4)
 */
static double compute_velocity_amplitude(int cell_idx)
{
    if (cell_idx < 1 || cell_idx > M) return 1e-4;

    double vmin = velocity_min[cell_idx];
    double vmax = velocity_max[cell_idx];
    double amp = 0.0;

    if (isfinite(vmin) && isfinite(vmax)) {
        amp = fmax(fmax(fabs(vmin), fabs(vmax)), 0.5 * (vmax - vmin));
    }
    if (!(amp > 0.0)) amp = fabs(velocity[cell_idx]);
    if (!(amp > 0.0)) amp = 1e-4;
    
    return amp;
}

/**
 * @brief Lateral mixing length Lm (Rutherford 1994, Fischer 1979)
 * @details Lm = B × [2 + 0.6 × √(Q_main/Q_trib)], clamped to [2B, 5B]
 */
static double estimate_lateral_mixing_length(double W, double vel_amp, double Q_trib, double A)
{
    if (!(W > 0.0)) return 0.0;
    
    double Q_main = fmax(fabs(vel_amp) * fmax(A, 0.0), 0.0);
    double q_ratio = (Q_trib > 1e-6 && Q_main > 0.0) ? sqrt(Q_main / Q_trib) : 0.0;
    double Lm = W * (2.0 + 0.6 * q_ratio);
    
    return fmin(fmax(Lm, 2.0 * W), 5.0 * W);
}

/**
 * @brief Dispersion coefficients using Fortran-heritage recursive formula
 * @details Savenije (2012) Eq. 5.72 Van der Burgh dispersion model:
 *   
 *   The steady-state dispersion equation (Savenije 2012, §5.4):
 *   dD/dx = -K × Q_f / A   where Q_f is freshwater discharge (positive downstream)
 *   
 *   Discrete form moving upstream (i increases towards river):
 *   D[i] = D[i-1] × (1 - beta × (exp(Δx/LCD) - 1))
 *   where beta = K × LCD × |Q_f| / (D[i-1] × A[i-1])
 *   
 *   Physical interpretation: Freshwater discharge REDUCES dispersion upstream
 *   because the salinity gradient weakens. Higher |Q_f| → faster decay → less salt intrusion.
 */
void Dispcoef(int t)
{
    ensure_mouth_geometry();

    double D0 = compute_mouth_dispersion(t);
    
    // Update scale-dependent dispersion bounds based on D0 (Reviewer feedback)
    // This ensures bounds adapt to estuary scale (Gisen et al. 2015)
    DISPERSION_MIN_DYNAMIC = fmax(0.1 * D0, dispersion_abs_min);  // 10% of D0
    DISPERSION_MAX_DYNAMIC = fmin(5.0 * D0, dispersion_abs_max);  // 5× D0
    
    // Total Q_fresh for logging only (at mouth)
    double Q_fresh_total = aggregate_freshwater_discharge(t);

    // --- Discharge-dependent K correction (Savenije 2005, 2012) ---
    // K_eff = K_geometric × (N / N_ref)^K_DISCHARGE_SENSITIVITY
    // Wet season (high Q → high N): K increases → faster salt decay
    // Dry season (low Q → low N):   K decreases → more salt intrusion
    double K_discharge_factor = 1.0;
    if (K_DISCHARGE_SENSITIVITY > 0.0) {
        double N_now = compute_canter_cremers_number(t);
        // Running exponential mean of N as reference (τ ≈ 365 days)
        static double N_ref_ema = -1.0;
        static long   N_step_count = 0;
        if (N_ref_ema < 0.0) {
            N_ref_ema = N_now;
        } else {
            // ~87600 steps/year at dt=360s (240 steps/day × 365)
            double ema_alpha = 1.0 / 87600.0;
            N_ref_ema = N_ref_ema * (1.0 - ema_alpha) + N_now * ema_alpha;
        }
        N_step_count++;
        // Apply only after warmup (~100 days = 24000 steps)
        if (N_step_count > 24000 && N_ref_ema > 0.0) {
            K_discharge_factor = pow(N_now / N_ref_ema, K_DISCHARGE_SENSITIVITY);
            // Safety clamp: avoid extreme values
            K_discharge_factor = fmin(fmax(K_discharge_factor, 0.2), 5.0);
        }
    }
    
    const int transition_center = index_2 > 0 ? index_2 : (M / 2);
    const int transition_width = segment_transition_width;

    disp[0] = disp[1] = D0;

    for (int i = 2; i <= M; ++i) {
        // Smooth transition for LC consistent with init.c (MAINT-001A fix)
        double segment_weight = 0.5 * (1.0 + tanh((double)(i - transition_center) / fmax((double)transition_width, 1.0)));
        double LC = (1.0 - segment_weight) * LC1 + segment_weight * LC2;
        if (!(LC > 0.0)) LC = 1.0;

        double K = compute_van_der_burgh(LC) * K_discharge_factor;
        // Clamp K after discharge correction to physical bounds [0, 1]
        K = fmin(fmax(K, 0.0), 1.0);
        double LCD = compute_local_convergence_length(i, LC);
        
        // Use numerical epsilon floor only (no non-physical mouth-area flooring).
        double A_prev = fmax(totalArea[i - 1], 1e-12);
        
        // Use physically meaningful freshwater magnitude, mode-independent.
        double Q_local = dispersion_freshwater_magnitude(t, i);
        
        // Savenije (2012) Eq. 5.72: beta MUST be positive for decay
        // beta = K × LCD × |Q_f| / (D × A)
        // When Q_local > 0, dispersion DECREASES moving upstream (proper salt flushing)
        double beta = (disp[i-1] > 0.0 && A_prev > 0.0) 
                    ? (K * LCD * Q_local) / (disp[i-1] * A_prev) 
                    : 0.0;
        
        double stretch = exp((double)DELXI / fmax(LCD, 1.0));
        
        // D[i] = D[i-1] × (1 - beta × (stretch - 1))
        // With positive beta, this ensures D[i] < D[i-1] (dispersion decays upstream)
        double delta = disp[i-1] * (1.0 - beta * (stretch - 1.0));

        // Enforce physical bounds on dispersion
        // Gisen (2015): D should be bounded by estuary scale. But the global
        // floor (0.1×D₀) is only appropriate in the tidal zone. Above the tidal
        // limit, the floor should decay with distance from the mouth, reflecting
        // the transition from tidal to fluvial dispersion (Fischer 1979).
        double dist_m = (double)(i * DELXI);
        double D_min_local = fmax(
            0.1 * D0 * exp(-dist_m / fmax(LC, 1.0)),
            dispersion_abs_min
        );
        
        if (delta < D_min_local) {
            delta = D_min_local;
        } else if (delta > DISPERSION_MAX_DYNAMIC) {
            delta = DISPERSION_MAX_DYNAMIC;
        }
        
        disp[i] = delta;
    }

    static int first = 1;
    if (first && debug_level > 0) {
        PRINTF_INIT("🌊 Dispersion: B₀=%.0fm h₀=%.1fm D₀=%.0fm²/s K₁=%.2f K₂=%.2f Q_f(total)=%.1fm³/s\n",
                    mouth_width_ref, mouth_depth_ref, D0,
                    compute_van_der_burgh(LC1), compute_van_der_burgh(LC2), Q_fresh_total);
        first = 0;
    }
}

/**
 * @brief Apply Boundary Conditions
 */
void Openbound(double* co, int s, int t)
{
    if (M < 3) return;

    const double dx = (double)DELXI;
    const double dt = (double)DELTI;
    const double clb = (downstreamBC[s].data && downstreamBC[s].dataSize > 0)
        ? interpolateInputData(t,
                               downstreamBC[s].time,
                               downstreamBC[s].data,
                               downstreamBC[s].dataSize,
                               false,
                               variableNames[s],
                               v[s].clb)
        : v[s].clb;

    double cub_val = (upstreamBC[s].data && upstreamBC[s].dataSize > 0)
        ? interpolateInputData(t,
                               upstreamBC[s].time,
                               upstreamBC[s].data,
                               upstreamBC[s].dataSize,
                               false,
                               variableNames[s],
                               v[s].cub)
        : v[s].cub;

    // Apply scaling for SPM upstream if override is active
    if (s == SPM && spm_upstream_scale != 1.0 && spm_upstream_scale > 0.0) {
        cub_val *= spm_upstream_scale;
    }
    const double cub = cub_val;

    // =========================================================================

    // --- Downstream Boundary (Mouth, cell j=1) ---
    {
        const double co1_prev = co[1];
        const double co3_prev = co[3];
        const double u_mouth = isfinite(velocity[2]) ? velocity[2] : 0.0;
        const double A_old = (isfinite(totalAreaOld[1]) && totalAreaOld[1] > 0.0) ? totalAreaOld[1] : totalArea[1];
        const double vol_old = fmax(A_old * dx, 1e-6);
        const double vol_new = fmax(totalArea[1] * dx, 1e-6);

        // =====================================================================
        // CFL-bounded upwind advection BC (parameter-free, all tracers)
        // =====================================================================
        // Standard 1st-order finite-volume upwind scheme applied at the open
        // boundary.  The blending coefficient α = min(|u|·dt/dx, 1) ensures
        // CFL stability and smooth transition through slack water (α→0).
        //
        //   Flood (u ≥ 0): upwind source = ocean boundary (clb)
        //                   co[1] = (1 - α)·co_prev + α·clb
        //   Ebb   (u < 0): upwind source = interior cell 3
        //                   co[1] = (1 - α)·co_prev + α·co3
        //
        // Properties:
        //   • Strong flood → mouth ≈ clb  (α ≈ 0.09/step → 99.7% per cycle)
        //   • Slack water → no change     (α ≈ 0, prevents ratchet/salt plug)
        //   • CFL-bounded → unconditionally stable
        // =====================================================================
        if (u_mouth >= 0.0) {
            const double alpha = fmin(u_mouth * dt / dx, 1.0);
            co[1] = co1_prev + alpha * (clb - co1_prev);
        } else {
            const double alpha = fmin(-u_mouth * dt / dx, 1.0);
            co[1] = co1_prev + alpha * (co3_prev - co1_prev);
        }

        // Diffusion/Relaxation (Flood only): Nudge towards ocean value (clb)
        // Characteristic timescale: tau_diff ~= dx^2 / D_mouth
        // Characteristic relaxation: beta = min(dt/tau_diff, 1.0)
        // 
        // IMPORTANT: Ebb (u < 0) uses PURE upwind advection to allow mass to exit
        // freely without being "plugged" by clb relaxation.
        if (u_mouth >= 0.0) {
            last_boundary_flux_mouth_disp_in[s] = 0.0;
            last_boundary_flux_mouth_disp_out[s] = 0.0;
            const double co_after_adv = co[1];
            const double D_mouth = (isfinite(disp[1]) && disp[1] > 0.0) ? disp[1] : 0.0;
            if (D_mouth > 0.0) {
                const double tau_diff = (dx * dx) / D_mouth;
                const double beta = fmin(fmax(dt / fmax(tau_diff, dt), 0.0), 1.0);
                co[1] = co_after_adv + beta * (clb - co_after_adv);

                const double disp_rate = (co[1] - co_after_adv) * vol_new / dt;
                last_boundary_flux_mouth_disp_in[s] = (disp_rate > 0.0) ? disp_rate : 0.0;
                last_boundary_flux_mouth_disp_out[s] = (disp_rate < 0.0) ? -disp_rate : 0.0;
            }
        } else {
            // Ebb: No dispersion/relaxation at boundary to prevent salt plug
            last_boundary_flux_mouth_disp_in[s] = 0.0;
            last_boundary_flux_mouth_disp_out[s] = 0.0;
        }

        {
            const double adv_rate = (co[1] - co1_prev) * vol_new / dt;
            last_boundary_flux_mouth_adv_in[s] = (adv_rate > 0.0) ? adv_rate : 0.0;
            last_boundary_flux_mouth_adv_out[s] = (adv_rate < 0.0) ? -adv_rate : 0.0;
            if (last_boundary_flux_mouth_disp_in[s] < 0.0 || !isfinite(last_boundary_flux_mouth_disp_in[s])) {
                last_boundary_flux_mouth_disp_in[s] = 0.0;
            }
            if (last_boundary_flux_mouth_disp_out[s] < 0.0 || !isfinite(last_boundary_flux_mouth_disp_out[s])) {
                last_boundary_flux_mouth_disp_out[s] = 0.0;
            }
        }

        if (!isfinite(co[1])) {
            co[1] = fmax(clb, 0.0);
        } else if (co[1] < 0.0) {
            co[1] = 0.0;
        }

        if (vol_new > 0.0) {
            last_boundary_flux_mouth[s] = (co[1] * vol_new - vol_old * co1_prev) / dt;
        }
    }

    // --- Upstream Boundary (River, upstream odd tracer cell) ---
    {
        const int i_up = upstream_tracer_index();
        const int i_in = (i_up >= 3) ? (i_up - 2) : i_up;
        const int face_interior = (i_up >= 2) ? (i_up - 1) : 2;
        const int face_boundary = ((M % 2) != 0) ? face_interior : fmin(i_up + 1, M);
        double u_up = velocity[face_boundary];
        if (!isfinite(u_up)) {
            u_up = velocity[face_interior];
        }
        if (!isfinite(u_up)) {
            u_up = 0.0;
        }
        const double coup_prev = co[i_up];
        const double Aup_old = (isfinite(totalAreaOld[i_up]) && totalAreaOld[i_up] > 0.0) ? totalAreaOld[i_up] : totalArea[i_up];
        const double vup_old = fmax(Aup_old * dx, 1e-6);
        const double vup_new = fmax(totalArea[i_up] * dx, 1e-6);

        if (u_up >= 0.0) {
            co[i_up] = co[i_up] - (co[i_up] - co[i_in]) * u_up * dt / dx;
        } else {
            co[i_up] = co[i_up] - (cub - co[i_up]) * u_up * dt / dx;
        }

        if (!isfinite(co[i_up])) {
            co[i_up] = fmax(cub, 0.0);
        } else if (co[i_up] < 0.0) {
            co[i_up] = 0.0;
        }

        // Keep even-index ghost consistent for output compatibility when M is even.
        if (i_up != M) {
            co[M] = co[i_up];
        }

        if (vup_new > 0.0) {
            last_boundary_flux_upstream[s] = (co[i_up] * vup_new - vup_old * coup_prev) / dt;
        }
    }
}

/**
 * @brief Tributary mass injection (Fischer 1979 / Rutherford 1994)
 * 
 * Two modes available via `tributary_mass_injection_mode`:
 *   - FISCHER: Single-cell exponential relaxation mixing
 *   - RUTHERFORD: Multi-cell lateral spread over [2B, 5B] mixing length
 * 
 * Core physics: C_new = C_trib + (C_old - C_trib) × exp(-Q_trib/V_mix × dt)
 */
void applyTributarySourceTerms(double* co, int s, int t)
{
    if (!tributaryEnabled || numTributaries <= 0) return;
    // Skip auto-calculated variables.
    if (s == pCO2 || s == PH || s == CO2) return;

    // Ablation toggle: skip lateral phytoplankton mixing from urban tributaries.
    // When active, the Fischer/Rutherford mixing step is skipped entirely for
    // Phy1/Phy2 at urban tributaries — no dilution, no seeding — so the river
    // concentration is unaffected, analogous to how pCO2/PH/CO2 are skipped.
    // Set CGEM_SKIP_LATERAL_PHYTO=1 to activate.
    static int skip_lateral_phyto = -1;
    if (skip_lateral_phyto < 0) {
        const char *env = getenv("CGEM_SKIP_LATERAL_PHYTO");
        skip_lateral_phyto = (env && atoi(env) > 0) ? 1 : 0;
        if (skip_lateral_phyto) {
            printf("[ABLATION] CGEM_SKIP_LATERAL_PHYTO=1: urban tributary Phy1/Phy2 mixing skipped (no dilution, no seeding)\n");
        }
    }

    const double dt = (double)DELTI;

    if (tributary_mass_injection_mode == TRIB_INJECTION_FISCHER) {
        // FISCHER (1979): Single-cell injection with exponential relaxation
        for (int trib = 0; trib < numTributaries; ++trib) {
            const int raw_cell = tributaries[trib].cellIndex;
            const int cell = tributary_continuity_index(raw_cell);
            if (cell < 1 || cell > M || cell >= MAXM) {
                continue;
            }

            const double Q_trib = calculatePhysicsBasedTributaryDischarge(t, trib);
            if (!(Q_trib > 0.0)) {
                continue;
            }

            // Get tributary concentration from time-series or static fallback
            double C_trib = 0.0;
            if (s == Sal) {
                C_trib = 0.0;
            } else if (tributaries[trib].chemicalData[s].timeArray &&
                       tributaries[trib].chemicalData[s].dataArray &&
                       tributaries[trib].chemicalData[s].dataSize > 0) {
                // Interpolate from time-series CSV data
                C_trib = interpolateInputData(
                    t,
                    tributaries[trib].chemicalData[s].timeArray,
                    tributaries[trib].chemicalData[s].dataArray,
                    tributaries[trib].chemicalData[s].dataSize,
                    false,
                    "Tributary concentration",
                    0.0);
            } else {
                // Fallback to static concentration
                C_trib = tributaries[trib].concentration[s];
            }

            // Skip lateral phytoplankton mixing at urban tributaries
            // equation applied, river concentration stays unchanged.
            if (skip_lateral_phyto && (s == Phy1 || s == Phy2) && tributaries[trib].is_urban) {
                continue;
            }

            // Validate concentration
            if (!isfinite(C_trib)) {
                diagnostics_fatal_tracer_error("applyTributarySourceTerms:C_trib", t, s, cell, C_trib, "nonfinite");
                return;
            }
            if (C_trib < 0.0) {
                diagnostics_fatal_tracer_error("applyTributarySourceTerms:C_trib", t, s, cell, C_trib, "negative");
                return;
            }
            if (s != Sal && C_trib == 0.0) {
                continue;
            }

            // Cell geometry
            const double A_cell = totalArea[cell];
            
            if (!(A_cell > 0.0)) {
                continue;
            }
            
            // Fischer (1979) Ch. 5: MASS SOURCE TERM for tributaries
            // Use the standard relaxation form (stable, bounded):
            //   dC/dt = (Q_trib / V_mix) * (C_trib - C)
            // Solution over timestep dt:
            //   C_new = C_trib + (C_old - C_trib) * exp(-(Q_trib/V_mix) * dt)
            // This automatically accounts for dilution via volume replacement and avoids
            // overshoot (alpha never exceeds 1.0).
            
            // Staggered-grid consistency:
            // Transport injects tributary mass into the same continuity-point control volume
            // as hydrodynamics (odd indices). That control volume spans two half-cells,
            // i.e., length ≈ 2*DELXI, except at the domain boundaries where Openbound()
            // uses half-cell volumes (DELXI).
            double dx_control = 2.0 * (double)DELXI;
            if (cell == 1 || cell == upstream_tracer_index()) {
                // Boundary control volume treated as half-cell.
                dx_control = (double)DELXI;
            }
            const double V_cell = A_cell * dx_control;  // Cell volume [m³]
            
            if (!(V_cell > 0.0)) {
                continue;
            }
            
            const double C_local = co[cell];
            if (!isfinite(C_local) || C_local < 0.0) {
                diagnostics_fatal_tracer_error("applyTributarySourceTerms:C_local", t, s, cell, C_local, (!isfinite(C_local) ? "nonfinite" : "negative"));
                return;
            }
            
            // Exponential relaxation mixing (Fischer 1979)
            const double mixing_rate = Q_trib / V_cell; // [1/s]
            const double alpha = 1.0 - exp(-mixing_rate * dt);
            double C_new = C_local + alpha * (C_trib - C_local);
            
            // Validate result (no silent clamping)
            if (!isfinite(C_new) || C_new < 0.0) {
                diagnostics_fatal_tracer_error("applyTributarySourceTerms:C_new", t, s, cell, C_new, (!isfinite(C_new) ? "nonfinite" : "negative"));
                return;
            }
            
            co[cell] = C_new;
            
            // Track mass flux for diagnostics [mass/s]
            if (s >= 0 && s < CHEM_COUNT) {
                const double mass_flux = Q_trib * C_trib;  // Mass added per second
                trib_mass_rate_last_step[s] += mass_flux;
            }
        }
        
        return;
    }

    // RUTHERFORD (1994): Multi-cell lateral spread
    // Modified spread for more distinct behavior from FISCHER mode
    typedef struct {
        int index;
        double weight;
        double area;
        double cell_length;
        double cross_fraction;
        double original_conc;
    } TributaryMixCell;

    // Exponential spread weights based on Rutherford (1994) lateral mixing theory
    // Weight decay: w[i] = exp(-i * dx / Lm) where Lm is mixing length
    // Pre-computed for efficiency assuming typical Lm ≈ 2-3 cells
    // These weights are normalized later to sum to 1.0
    // NOTE: Must remain a compile-time constant (MSVC C does not support VLAs).
    #ifndef SPREAD_COUNT
    #define SPREAD_COUNT 5
    #endif
    static const int spread_offsets[SPREAD_COUNT] = {0, 1, 2, 3, 4};
    // Exponential decay: exp(-0), exp(-0.5), exp(-1), exp(-1.5), exp(-2)
    // = 1.0, 0.607, 0.368, 0.223, 0.135 (normalized sum ≈ 2.33)
    static const double spread_weights_base[SPREAD_COUNT] = {1.0, 0.607, 0.368, 0.223, 0.135};



    for (int trib = 0; trib < numTributaries; ++trib) {
        const int raw_cell = tributaries[trib].cellIndex;
        const int cell = tributary_continuity_index(raw_cell);
        if (cell < 1 || cell > M || cell >= MAXM) {
            continue;
        }

        const double Q_trib = calculatePhysicsBasedTributaryDischarge(t, trib);
        if (!(Q_trib > 0.0)) {
            continue;
        }

        double C_trib = 0.0;
        if (s == Sal) {
            C_trib = 0.0;
        } else if (tributaries[trib].chemicalData[s].timeArray &&
                   tributaries[trib].chemicalData[s].dataArray &&
                   tributaries[trib].chemicalData[s].dataSize > 0) {
            C_trib = interpolateInputData(
                t,
                tributaries[trib].chemicalData[s].timeArray,
                tributaries[trib].chemicalData[s].dataArray,
                tributaries[trib].chemicalData[s].dataSize,
                false,
                "Tributary concentration",
                0.0);
        } else {
            C_trib = tributaries[trib].concentration[s];
        }

        // Skip lateral phytoplankton mixing at urban tributaries (Rutherford path).
        if (skip_lateral_phyto && (s == Phy1 || s == Phy2) && tributaries[trib].is_urban) {
            continue;
        }

        if (!isfinite(C_trib)) {
            diagnostics_fatal_tracer_error("applyTributarySourceTerms:C_trib", t, s, cell, C_trib, "nonfinite");
            return;
        }
        if (C_trib < 0.0) {
            diagnostics_fatal_tracer_error("applyTributarySourceTerms:C_trib", t, s, cell, C_trib, "negative");
            return;
        }
        if (s != Sal && C_trib == 0.0) {
            continue;
        }

        const double area_local = (cell >= 1 && cell <= M) ? totalArea[cell] : totalArea[1];
        double vel_amp_local = compute_velocity_amplitude(cell);
        if (!(vel_amp_local > 1e-6)) {
            vel_amp_local = 0.05; // Minimal tidal stirring [m/s]
        }
        const double width_local = (cell >= 1 && cell <= M && width[cell] > 1.0) ? width[cell] : width[1];
        const double mixing_length = estimate_lateral_mixing_length(width_local,
                                                                    vel_amp_local,
                                                                    Q_trib,
                                                                    area_local);

        TributaryMixCell mix_cells[SPREAD_COUNT];
        int mix_cell_count = 0;
        double weight_sum = 0.0;

        for (int idx = 0; idx < SPREAD_COUNT; ++idx) {
            // Use exponential decay weights (Rutherford 1994)
            const double weight = spread_weights_base[idx];
            if (!(weight > 0.0)) {
                continue;
            }

            const int offset = spread_offsets[idx];
            // Longitudinal mixing zone should extend downstream of the confluence
            // (Rutherford 1994; Fischer 1979). In this model, cell index increases
            // upstream (km from mouth), so downstream is toward smaller indices.
            const int target = cell - offset;
            if (target < 1 || target > M) {
                continue;
            }

            const double reach_start = (double)offset * (double)DELXI;
            if (reach_start >= mixing_length) {
                continue;
            }

            const double area_target = totalArea[target];
            if (!(area_target > 0.0)) {
                continue;
            }

            double cell_length = mixing_length - reach_start;
            cell_length = fmin(cell_length, (double)DELXI);
            const double min_length = 0.25 * (double)DELXI;
            if (cell_length < min_length) {
                cell_length = min_length;
            }

            double Q_total_local = Discharge(t, target, NULL);
            double magnitude_total = fabs(Q_total_local);
            if (!(magnitude_total > 1e-6)) {
                magnitude_total = Q_trib;
            }

            double cross_fraction = 0.0;
            const double denom = Q_trib + magnitude_total;
            if (denom > 1e-6) {
                cross_fraction = Q_trib / denom;
            } else {
                cross_fraction = 0.5;
            }
            cross_fraction = fmin(0.5, fmax(0.05, cross_fraction));

            mix_cells[mix_cell_count].index = target;
            mix_cells[mix_cell_count].weight = weight;
            mix_cells[mix_cell_count].area = area_target;
            mix_cells[mix_cell_count].cell_length = cell_length;
            mix_cells[mix_cell_count].cross_fraction = cross_fraction;
            {
                const double c0 = co[target];
                if (!isfinite(c0) || c0 < 0.0) {
                    diagnostics_fatal_tracer_error("applyTributarySourceTerms:original_conc", t, s, target, c0, (!isfinite(c0) ? "nonfinite" : "negative"));
                    return;
                }
                mix_cells[mix_cell_count].original_conc = c0;
            }
            mix_cell_count++;
            weight_sum += weight;
        }

        if (mix_cell_count == 0 || !(weight_sum > 0.0)) {
            continue;
        }

        double mass_rate_accumulated = 0.0;

        for (int idx = 0; idx < mix_cell_count; ++idx) {
            const double normalized_weight = mix_cells[idx].weight / weight_sum;
            if (!(normalized_weight > 0.0)) {
                continue;
            }

            const double engaged_fraction = mix_cells[idx].cross_fraction * normalized_weight;
            const double volume_effective = mix_cells[idx].area * mix_cells[idx].cell_length * fmax(engaged_fraction, 0.0);
            if (!(volume_effective > 0.0)) {
                continue;
            }

            const double C_local = mix_cells[idx].original_conc;

            const double Q_share = Q_trib * normalized_weight;
            const double mixing_rate = (volume_effective > 0.0)
                ? (Q_share / volume_effective)
                : 0.0;
            const double decay = exp(-mixing_rate * dt);

            double C_new = C_trib + (C_local - C_trib) * decay;
            if (!isfinite(C_new) || C_new < 0.0) {
                diagnostics_fatal_tracer_error("applyTributarySourceTerms:C_new", t, s, mix_cells[idx].index, C_new, (!isfinite(C_new) ? "nonfinite" : "negative"));
                return;
            }

            const double gradient = C_trib - C_local;
            if (gradient > 0.0 && C_new > C_trib) {
                C_new = C_trib;
            } else if (gradient < 0.0 && C_new < C_trib) {
                C_new = C_trib;
            }

            co[mix_cells[idx].index] = C_new;

            const double mass_added = (C_new - C_local) * volume_effective;
            mass_rate_accumulated += mass_added / dt;
        }

        if (s >= 0 && s < CHEM_COUNT) {
            trib_mass_rate_last_step[s] += mass_rate_accumulated;
        }
    }
}

// -----------------------------------------------------------------------------
// Staggered-grid consistency helpers
// -----------------------------------------------------------------------------

/**
 * @brief Map a tributary cell index to the continuity-point index used by hydrodynamics.
 *
 * Hydrodynamics injects tributary discharge into the continuity equation at odd indices
 * (see src/hydrodynamics.c compute_lateral_inflow_profile()). To preserve strict
 * mass/volume consistency on the staggered grid, transport should inject tributary mass
 * into the same continuity-point control volume.
 */
static int tributary_continuity_index(int cell)
{
    int idx = cell;
    if ((idx % 2) == 0) {
        if (idx > 1) {
            idx -= 1;
        } else if (idx < M) {
            idx += 1;
        }
    }
    return idx;
}

/**
 * @brief Upstream tracer control-volume index on odd continuity points.
 */
static int upstream_tracer_index(void)
{
    return ((M % 2) != 0) ? M : (M - 1);
}

/**
 * @brief Advection Scheme (TVD)
 */
void TVD(double* co, int s)
{
    (void)s; // scheme is identical across tracers; species handled by caller

    const double dx = (double)DELXI;
    const double dt = (double)DELTI;
    const int i_up = upstream_tracer_index();

    // Use a local copy for the sub-stepping to ensure strict conservation
    static double cold[MAXM + 1];
    static double fl_accum[MAXM + 1];
    
    // Determine CFL for sub-stepping
    double max_v = 1e-6;
    for (int i = 2; i <= M; i += 2) {
        if (fabs(velocity[i]) > max_v) max_v = fabs(velocity[i]);
    }
    // Advection CFL: u * dt / dx < 1.0 (for 2dx cells, u*dt/2dx < 1.0)
    double cfl_max = max_v * dt / (2.0 * dx); 
    int nsub = (int)ceil(cfl_max / 0.8);
    if (nsub < 1) nsub = 1;
    if (nsub > 20) nsub = 20;

    tvd_calls_total++;
    if (transport_time_context >= WARMUP) {
        tvd_calls_post_warmup++;
    }
    if (nsub > tvd_max_substeps_observed) {
        tvd_max_substeps_observed = nsub;
    }
    if (nsub > 1) {
        tvd_substep_activation_events++;
        if (transport_time_context >= WARMUP) {
            tvd_substep_activation_events_post_warmup++;
        }
    }

    const double dt_sub = dt / (double)nsub;
    const double cfl_sub = max_v * dt_sub / (2.0 * dx);

    for (int i = 0; i <= MAXM; i++) fl_accum[i] = 0.0;

    for (int sub = 0; sub < nsub; sub++) {
        for (int i = 1; i <= M; i++) cold[i] = co[i];
        
        // --- 1. Compute Fluxes ---
        // Faces are even indices 2, 4, ..., M1.
        // fl[0] and fl[M+1] are reserved for Openbound.
        for (int i = 0; i <= MAXM; i++) fl[i] = 0.0;

        for (int face = 2; face <= i_up - 1; face += 2) {
            const double vx = velocity[face];
            const int j = face - 1; // Left cell
            const int jr = face + 1; // Right cell
            
            double philen = 0.0;
            double f = 0.0;
            double rg = 0.0;
            double phi = 0.0;

            if (vx > 0.0) { // Upstream flow (+x direction)
                // For face=2 (j=1, jr=3), the upstream cell is j=1. The cell to its left is j-2 (which is -1).
                // The gradient is calculated using cold[j-2], cold[j], cold[jr].
                // If j=1, we don't have cold[j-2].
                // The original TVD code used j!=1 for rg calculation.
                // Here, we use j >= 3 to ensure cold[j-2] is valid.
                if (jr + 1 <= M) { // Check if cold[jr+1] (i.e., cold[j+2]) is valid
                    f = cold[jr] - cold[j]; // f = cold[j+2] - cold[j]
                    if (fabs(f) > TVD_FLUX_THRESHOLD && j >= 3) { // Check if cold[j-2] is valid
                        rg = (cold[j] - cold[j - 2]) / f;
                        phi = (2.0 - cfl_sub) / 3.0 + (1.0 + cfl_sub) / 3.0 * rg;
                        philen = fmax(0.0, fmin(2.0, fmin(2.0 * rg, phi)));
                        if (philen > 0.0) {
                            limiter_face_activation_events++;
                            if (transport_time_context >= WARMUP) {
                                limiter_face_activation_events_post_warmup++;
                            }
                        }
                    }
                }
                double C_face = cold[j] + 0.5 * (1.0 - cfl_sub) * philen * f;
                fl[face] = vx * totalArea[face] * C_face;
            } else if (vx < 0.0) { // Downstream flow (-x direction)
                // For face=M1 (j=M1-1, jr=M1+1), the upstream cell is jr=M1+1. The cell to its right is jr+2.
                // If jr=M1+1, we don't have cold[jr+2].
                // The original TVD code used j!=M2 for rg calculation.
                // Here, we use jr <= M-2 to ensure cold[jr+2] is valid.
                if (j - 1 >= 1) { // Check if cold[j-1] (i.e., cold[j-2]) is valid
                    f = cold[j] - cold[jr]; // f = cold[j] - cold[j+2]
                    if (fabs(f) > TVD_FLUX_THRESHOLD && jr <= M - 2) { // Check if cold[jr+2] is valid
                        rg = (cold[jr] - cold[jr + 2]) / f;
                        phi = (2.0 - cfl_sub) / 3.0 + (1.0 + cfl_sub) / 3.0 * rg;
                        philen = fmax(0.0, fmin(2.0, fmin(2.0 * rg, phi)));
                        if (philen > 0.0) {
                            limiter_face_activation_events++;
                            if (transport_time_context >= WARMUP) {
                                limiter_face_activation_events_post_warmup++;
                            }
                        }
                    }
                }
                double C_face = cold[jr] + 0.5 * (1.0 - cfl_sub) * philen * f;
                fl[face] = vx * totalArea[face] * C_face;
            }
            fl_accum[face] += fl[face];
        }

        // --- 2. Update Concentrations (proper FVM volume accounting) ---
        // Fortran-ported: TVD updates INTERIOR cells only (j=3..M-2).
        // Boundary cells (j=1, j=M1) are handled exclusively by Openbound,
        // which runs BEFORE TVD. This avoids double-counting at boundaries.
        // NOTE: Must account for volume change between old and new timestep.
        for (int j = 3; j <= i_up - 2; j += 2) {
            const double cell_len = 2.0 * dx;  // interior cells span 2*dx
            const double A_old_j = (isfinite(totalAreaOld[j]) && totalAreaOld[j] > 0.0) ? totalAreaOld[j] : totalArea[j];
            const double Vol_old = fmax(A_old_j * cell_len, 1e-6);
            const double Vol_new = fmax(totalArea[j] * cell_len, 1e-6);
            const double flux_div = fl[j + 1] - fl[j - 1];
            co[j] = (cold[j] * Vol_old - dt_sub * flux_div) / Vol_new;
        }
    }

    // Average fluxes for diagnostis
    for (int i = 0; i <= MAXM; i++) fl[i] = fl_accum[i] / (double)nsub;
}

/**
 * @brief Implicit Crank-Nicolson Dispersion Scheme (ported from Fortran)
 *
 * Replaces the previous explicit diffusion with an implicit tridiagonal solver.
 * This is unconditionally stable regardless of dt, D, or dx — no CFL limit.
 *
 * Fortran reference: CGEM_Transport.f90 SUBROUTINE Calculate_Dispersion (L370-460)
 *
 * Boundary conditions: Dirichlet (a=0, b=1, c=0, d=conc) at j=1 and j=M1.
 * Interior: Crank-Nicolson coefficients with tridiag solve.
 */
void Disp(double* co, int s)
{
    (void)s; // diffusion operator independent of tracer identity

    const double dt = (double)DELTI;
    const double dx = (double)DELXI;
    const int i_up = upstream_tracer_index();
    const int i_up_inner = (i_up >= 3) ? (i_up - 2) : 1;

    // Tridiag arrays (only odd indices used: 1, 3, 5, ..., M1)
    // We index by the actual grid index for clarity.
    static double ta[MAXM + 1];   // sub-diagonal
    static double tb[MAXM + 1];   // diagonal
    static double tc[MAXM + 1];   // super-diagonal
    static double td[MAXM + 1];   // RHS
    static double tgam[MAXM + 1]; // scratch for Thomas algorithm
    static double cold[MAXM + 1]; // old concentrations

    // Store old values
    for (int j = 1; j <= i_up; j += 2) {
        cold[j] = co[j];
    }

    // --- Boundary conditions (Dirichlet: hold current concentration) ---
    // Fortran: a(1)=0, b(1)=1, c(1)=0, di(1)=conc(1)
    ta[1]  = 0.0;
    tb[1]  = 1.0;
    tc[1]  = 0.0;
    td[1]  = co[1];

    ta[i_up] = 0.0;
    tb[i_up] = 1.0;
    tc[i_up] = 0.0;
    td[i_up] = co[i_up];

    // --- Interior Crank-Nicolson coefficients ---
    // Fortran: alpha = DELTI / (8 * DELXI * DELXI)
    const double cn_alpha = dt / (8.0 * dx * dx);

    for (int i = 3; i <= i_up_inner; i += 2) {
        // g1 = D(i-1) * A(i-1) / A(i)   (left face contribution)
        // g2 = D(i+1) * A(i+1) / A(i)   (right face contribution)
        const double D_left  = (isfinite(disp[i-1]) && disp[i-1] > 0.0) ? disp[i-1] : 0.0;
        const double D_right = (isfinite(disp[i+1]) && disp[i+1] > 0.0) ? disp[i+1] : 0.0;
        const double A_i = fmax(totalArea[i], 1e-12);
        const double g1 = D_left  * fmax(totalArea[i-1], 1e-12) / A_i;
        const double g2 = D_right * fmax(totalArea[i+1], 1e-12) / A_i;

        ta[i] = -g1 * cn_alpha;
        tc[i] = -g2 * cn_alpha;
        tb[i] = 1.0 + (g1 + g2) * cn_alpha;
        double r_coeff = 1.0 - (g1 + g2) * cn_alpha;

        // Fortran: di(i) = -c(i)*conc(i+2) + r(i)*conc(i) - a(i)*conc(i-2)
        td[i] = -tc[i] * co[i+2] + r_coeff * co[i] - ta[i] * co[i-2];
    }

    // --- Thomas algorithm (tridiag solver) ---
    // Fortran: forward sweep then back-substitution
    // Note: the system is indexed by odd nodes: 1, 3, 5, ..., M1
    // In the Thomas algorithm, "previous" means j-2 (odd spacing)
    double bet = tb[1];
    if (fabs(bet) < 1e-30) {
        printf("Error Disp: b[1]=0.0\n");
        return;
    }
    co[1] = td[1] / bet;

    // Forward sweep
    for (int j = 3; j <= i_up; j += 2) {
        tgam[j] = tc[j-2] / bet;
        bet = tb[j] - ta[j] * tgam[j];
        if (fabs(bet) < 1e-30) {
            // Degenerate row — keep old value
            co[j] = cold[j];
            continue;
        }
        co[j] = (td[j] - ta[j] * co[j-2]) / bet;
    }

    // Back-substitution
    for (int j = i_up_inner; j >= 1; j -= 2) {
        co[j] = co[j] - tgam[j+2] * co[j+2];
    }

    // Keep even-index points consistent (interpolation for outputs)
    for (int i = 2; i <= i_up - 1; i += 2) {
        co[i] = 0.5 * (co[i - 1] + co[i + 1]);
    }

    // Upstream ghost
    if (i_up != M && i_up >= 1 && i_up <= M) {
        co[M] = co[i_up];
    }
}

void Boundflux(int s)
{
    for(int j = 2; j <= M1; j += 2) {
        v[s].advflux[j] += fl[j] * DELTI;
        v[s].concflux[j] += v[s].c[j];
    }
}

void Transport(int t)
{
    if (t == 0) {
        tvd_calls_total = 0;
        tvd_calls_post_warmup = 0;
        tvd_substep_activation_events = 0;
        tvd_substep_activation_events_post_warmup = 0;
        tvd_max_substeps_observed = 1;
        limiter_face_activation_events = 0;
        limiter_face_activation_events_post_warmup = 0;
        diffusion_substep_activation_events = 0;
        diffusion_substep_activation_events_post_warmup = 0;
    }
    transport_time_context = t;

    // IMPORTANT: Dispersion depends on instantaneous geometry/forcing.
    // Update every timestep to avoid phase lag between hydrodynamics and transport.
    Dispcoef(t);

    // Define as a nested macro to keep MSVC C compatibility (no nested functions).
    #define VALIDATE_SPECIES_FIELD(STAGE, TLOCAL, FIELD, SLOCAL) do { \
        for (int _i = 1; _i <= M; ++_i) { \
            const double _val = (FIELD)[_i]; \
            if (!isfinite(_val) || _val < 0.0) { \
                diagnostics_fatal_tracer_error((STAGE), (TLOCAL), (SLOCAL), _i, _val, (!isfinite(_val) ? "nonfinite" : "negative")); \
                return; \
            } \
        } \
    } while (0)

    for (int s = 0; s < MAXV; s++) {
        if (v[s].env != 1) continue;
        if (s == pCO2 || s == PH || s == CO2) continue;

        VALIDATE_SPECIES_FIELD("Transport:entry", t, v[s].c, s);

        applyTributarySourceTerms(v[s].c, s, t);
        VALIDATE_SPECIES_FIELD("Transport:applyTributarySourceTerms", t, v[s].c, s);

        // Fortran call order: Openbound → TVD → Disp
        // Openbound handles boundary cells (j=1, M1) FIRST.
        // TVD then updates interior cells only (j=3..M-2), skipping boundaries.
        // Disp (implicit Crank-Nicolson) solves for all cells with Dirichlet BCs.
        Openbound(v[s].c, s, t);
        VALIDATE_SPECIES_FIELD("Transport:Openbound", t, v[s].c, s);
        TVD(v[s].c, s);
        VALIDATE_SPECIES_FIELD("Transport:TVD", t, v[s].c, s);
        Disp(v[s].c, s);
        VALIDATE_SPECIES_FIELD("Transport:Disp", t, v[s].c, s);

        if (t >= WARMUP) {
            for(int i = 1; i <= M; i++) v[s].avg[i] += v[s].c[i];
            Boundflux(s);
        }
    }

    #undef VALIDATE_SPECIES_FIELD
}

void export_limiter_events_summary(const char *path)
{
    if (!path || path[0] == '\0') {
        return;
    }

    ensure_directory_for_path(path);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        if (debug_level >= DEBUG_LEVEL_WARNING) {
            printf("⚠️  Could not write limiter diagnostics CSV: %s\n", path);
        }
        return;
    }

    const double tvd_substep_activation_fraction_post_warmup =
        (tvd_calls_post_warmup > 0)
            ? ((double)tvd_substep_activation_events_post_warmup / (double)tvd_calls_post_warmup)
            : 0.0;
    const double limiter_face_activation_rate_per_tvd_call =
        (tvd_calls_post_warmup > 0)
            ? ((double)limiter_face_activation_events_post_warmup / (double)tvd_calls_post_warmup)
            : 0.0;

    fprintf(fp,
            "tvd_calls_total_count,tvd_calls_post_warmup_count,tvd_substep_activation_events_count,tvd_substep_activation_events_post_warmup_count,tvd_substep_activation_fraction_post_warmup,tvd_max_substeps_observed_count,limiter_face_activation_events_count,limiter_face_activation_events_post_warmup_count,limiter_face_activation_rate_per_tvd_call_count_per_call,diffusion_substep_activation_events_count,diffusion_substep_activation_events_post_warmup_count,unit_notes\n");
    fprintf(fp,
            "%ld,%ld,%ld,%ld,%.10f,%d,%ld,%ld,%.10f,%ld,%ld,counts_and_dimensionless_fractions\n",
            tvd_calls_total,
            tvd_calls_post_warmup,
            tvd_substep_activation_events,
            tvd_substep_activation_events_post_warmup,
            tvd_substep_activation_fraction_post_warmup,
            tvd_max_substeps_observed,
            limiter_face_activation_events,
            limiter_face_activation_events_post_warmup,
            limiter_face_activation_rate_per_tvd_call,
            diffusion_substep_activation_events,
            diffusion_substep_activation_events_post_warmup);

    fclose(fp);
}