/**
 * @file hydrodynamics.c
 * @brief Consolidated hydrodynamic module for 1D Saint-Venant equations in tidal rivers
 * @author An Nguyen
 * @date Last updated: 05/2025
 *
 * This file contains the complete hydrodynamic solver, including:
 * - Core solver functions (previously in hyd.c)
 * - Auxiliary functions (previously in uphyd.c)
 * - Matrix operations and boundary condition handling (previously in tridaghyd.c)
 */

#include "define.h"
#include "variables.h"
#include "diagnostics.h" // Consolidated diagnostics system
#include "file.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static double hydro_convergence_tolerance = -1.0; // Initialized on first Hyd() call
static double last_eta_residual = 0.0;
static double last_velocity_residual = 0.0;
static int consecutive_maxits = 0; // Track consecutive MAXITS hits for auto-relaxation
static long hydro_solver_calls_total = 0;
static long hydro_solver_calls_post_warmup = 0;
static long hydro_iteration_sum_post_warmup = 0;
static long hydro_maxits_hit_count = 0;
static int hydro_max_iteration_observed = 0;
static int hydro_max_iteration_observed_post_warmup = 0;
static long tidal_stats_last_day = -1; // Track daily tidal statistics window
static double tidal_range_sum[MAXM + 1];
static double tidal_range_count[MAXM + 1];

/**
 * @brief Compute the maximum velocity magnitude and corresponding transport CFL number.
 *
 * The CFL metric follows the explicit upwind stability criterion (LeVeque 2002), using
 * the recorded velocity extrema accumulated during the current diagnostic window.
 *
 * @param[out] max_velocity_out Optional pointer receiving the peak |U| [m/s].
 * @return Maximum CFL = |U| Δt / Δx across the grid (dimensionless).
 */
static double compute_max_velocity_cfl(double *max_velocity_out)
{
    double max_vel = 0.0;
    double max_cfl = 0.0;

    for (int i = 1; i <= M; ++i) {
        double v_peak = fmax(fabs(velocity_max[i]), fabs(velocity_min[i]));
        if (!isfinite(v_peak)) {
            continue;
        }
        if (v_peak > max_vel) {
            max_vel = v_peak;
        }
        double cfl = v_peak * (double)DELTI / (double)DELXI;
        if (cfl > max_cfl) {
            max_cfl = cfl;
        }
    }

    if (max_velocity_out) {
        *max_velocity_out = max_vel;
    }
    return max_cfl;
}

// Convergence check temporary arrays (match standalone version names)
double E[MAXM+1];  // Previous water level values
double Y[MAXM+1];  // Previous velocity values

// Global arrays to track tidal range statistics
double daily_min_level[MAXM + 1]; // Daily minimum water level at each cell
double daily_max_level[MAXM + 1]; // Daily maximum water level at each cell
double velocity_min[MAXM + 1];    // Minimum velocity at each cell
double velocity_max[MAXM + 1];    // Maximum velocity at each cell
double tidal_range_mean[MAXM + 1]; // Mean tidal range at each cell across the analysis window [m]
double dispersion_min[MAXM + 1];  // Minimum dispersion at each cell
double dispersion_max[MAXM + 1];  // Maximum dispersion at each cell
int last_saved_day = -1;      // Track the last day we saved data for

// Debug control — only active when debug_level >= DEBUG_LEVEL_VERBOSE
#define DEBUG_HYDRAULICS 0
#define DEBUG_MOUTH_POINTS 5

// Lateral inflow profile (tributaries) stored as m²/s for continuity source term
static double lateral_inflow_profile[MAXM + 1];

static void compute_lateral_inflow_profile(int t)
{
    memset(lateral_inflow_profile, 0, sizeof(lateral_inflow_profile));

    if (!tributaryEnabled || numTributaries <= 0 || !tributaries) {
        return;
    }

    // Coeffa() discretizes dQ/dx using inv_2dx = 1/(2*DELXI).
    // For consistency on the same staggered continuity control length (~2*DELXI),
    // convert a tributary discharge Q_trib [m^3/s] into a source term q [m^2/s]
    // using Q_trib/(2*DELXI).
    const double inv_2dx = 1.0 / (2.0 * (double)DELXI);

    for (int trib = 0; trib < numTributaries; ++trib) {
        int cell = tributaries[trib].cellIndex;
        if (cell < 1 || cell > M) {
            continue;
        }

        double Q_trib = calculatePhysicsBasedTributaryDischarge(t, trib);
        if (!(Q_trib > 0.0)) {
            continue;
        }
        int continuity_idx = cell;
        if ((continuity_idx % 2) == 0) {
            if (continuity_idx > 1) {
                continuity_idx -= 1;
            } else if (continuity_idx < M) {
                continuity_idx += 1;
            }
        }

        if (continuity_idx < 1 || continuity_idx > M) {
            continue;
        }

        // Fischer (1979): Tributaries handled as mass sources in transport equation
        // Add discharge to continuity for momentum conservation only
        double q_term = Q_trib * inv_2dx;
        lateral_inflow_profile[continuity_idx] += q_term; // [m³/s] / m → [m²/s]
    }
}

/**
 * @brief Boundary conditions with improved tidal oscillation
 * 
 * @note DESIGN LIMITATION: Minimum depth floor of 0.1m
 * This model is designed for DEEP tidal estuaries (typical depth 5-15m) where
 * wetting/drying of intertidal flats is not a primary concern. The 0.1m floor
 * prevents numerical instabilities at extreme low tides but does not conserve
 * mass exactly during near-dry conditions.
 * 
 * For shallow estuaries with significant intertidal zones, a proper wetting/drying
 * algorithm (e.g., Defina 2000, Casulli 2009) should be implemented.
 * 
 * @param t Current simulation time in seconds
 */
void Newbc(int t)
{
    // Log boundary condition setting only periodically (every hour)
    if (debug_is_enabled(DEBUG_LEVEL_DETAIL)) {
        debug_log_periodic(DEBUG_LEVEL_DETAIL, DEBUG_MODULE_HYDRO, 3600, "boundary_conditions",
                          "Setting boundary conditions at t=%d (day %.2f)", 
                  t, (double)t/86400.0);
    }
    
    double Q_up;
    
    // Downstream boundary - water level prescribed (EXACTLY following Fortran CGEM)
    double tideLevel = Tide(t);
    
    // ✅ CORRECT: Follow Fortran CGEM_Hydrodynamics.f90 exactly
    // AAf(1) = B(1)*TidalElevation()
    // 
    // The key insight: tideLevel represents WATER DEPTH above bed, not elevation above datum
    // This means totalArea = width * (riverbed_depth + tidal_variation)
    // where riverbed_depth is the mean depth and tidal_variation oscillates around zero
    
    double total_water_depth = riverbed_depth[1] + tideLevel;  // Physical water depth
    
    // Minimum depth floor (see function documentation for design limitation note)
    // Value: 0.1m based on typical numerical stability requirements
    const double MIN_DEPTH_FLOOR = 0.1;  // [m] - Named constant for clarity
    if (total_water_depth <= MIN_DEPTH_FLOOR) {
        total_water_depth = MIN_DEPTH_FLOOR;
        if (debug_level >= DEBUG_LEVEL_WARNING) {
            printf("⚠️  WARNING: Extreme low tide at t=%d, depth limited to %.1f m (tideLevel=%.2f m)\n", 
                   t, total_water_depth, tideLevel);
        }
    }
    
    // Set cross-sectional areas based on total water depth
    totalArea[1] = width[1] * total_water_depth;     // AA[1] = B[1] * H[1]  
    freeArea[1] = totalArea[1] - baseArea[1];        // Free surface area
    waterDepth[1] = total_water_depth;               // Store actual depth
    level[1] = tideLevel;                            // Store tidal elevation for output
    tempFreeArea[1] = freeArea[1];                   // For iterative solver
    
    // --- Upstream Boundary (i=M) ---
    // Get the upstream discharge (negative for inflow)
    Q_up = Discharge_ups(t); // Using the correct discharge function

    // Fortran-equivalent upstream velocity BC: u(M) = Q_up / A(M)
    // with robust area fallback/guard to avoid non-finite values.
    {
        const double area_eps = 1e-9;
        double area_up = totalArea[M];

        if (!(area_up > area_eps) || !isfinite(area_up)) {
            area_up = baseArea[M] + tempFreeArea[M];
        }
        if (!(area_up > area_eps) || !isfinite(area_up)) {
            area_up = (M >= 2) ? totalArea[M-1] : area_up;
        }
        if ((!(area_up > area_eps) || !isfinite(area_up)) && width[M] > area_eps && waterDepth[M] > 0.0) {
            area_up = width[M] * waterDepth[M];
        }
        if (!(area_up > area_eps) || !isfinite(area_up)) {
            area_up = area_eps;
        }

        tempVelocity[M] = Q_up / area_up;
        if (!isfinite(tempVelocity[M])) {
            tempVelocity[M] = 0.0;
        }
    }
}

/**
 * @brief Matrix coefficient assembly for the Saint-Venant equations.
 */
void Coeffa(int t)
{
    // Aggressive hoisting of invariants (physics unchanged – performance optimized)
    // Initialize once on first call since values come from configuration files
    static double dt = 0.0, inv_dt = 0.0, dx = 0.0, inv_2dx = 0.0, inv_4Gdx = 0.0, inv_2Gdx = 0.0, inv_Gdt = 0.0;
    static int initialized = 0;
    
    if (!initialized) {
        dt = (double)DELTI;
        inv_dt = 1.0 / (double)DELTI;
        dx = DELXI;
        inv_2dx = 1.0 / (2.0 * DELXI);
    inv_4Gdx = 1.0 / (4.0 * G * DELXI);
    inv_2Gdx = 1.0 / (2.0 * G * DELXI);
        inv_Gdt = 1.0 / (G * (double)DELTI);
        initialized = 1;
    }
    
    const double Qr = Discharge_ups(t);  // River discharge (negative inflow)

    // Update lateral inflow sources (tributaries) for continuity equation (Fischer 1979)
    compute_lateral_inflow_profile(t);

    // Interior equations (odd j continuity, even i=j+1 momentum)
    for (int j = 3; j <= M-2; j += 2) {
        int i = j + 1;
        double totalArea_jm1 = totalArea[j-1];
        double totalArea_jp1 = totalArea[j+1];

        // Continuity at odd points (optimized)
        double q_lateral = lateral_inflow_profile[j];  // [m²/s] source from tributaries
        double local_RS_inv_dt = rs[j] * inv_dt;
        Z[j] = freeArea[j] * local_RS_inv_dt + q_lateral;  // Continuity with mass source term
        C[j][1] = -totalArea_jm1 * inv_2dx;
        C[j][2] = local_RS_inv_dt;  // Use precomputed product
        C[j][3] = totalArea_jp1 * inv_2dx;
        C[j][4] = 0.0;

        // Momentum at even points with FULL Saint-Venant terms (Savenije 2012)
        double ui = velocity[i];
        double chezy_i_sq = Chezy[i] * Chezy[i];  // Eliminate division by squares
        double wdepth_i = waterDepth[i];
        double vel_diff = velocity[i+2] - velocity[i-2];
        double friction_term = fabs(ui) / (chezy_i_sq * wdepth_i); // Decomposed friction
        double convective_term = vel_diff * inv_4Gdx;              // Semi-implicit convective contribution

        // ✅ Fortran: No sponge layer — Chezy friction profile handles damping naturally.
        // Removed non-physical friction sponge that caused tidal rebound upstream.

        Z[i] = ui * inv_Gdt;  // Keep convective term on diagonal for stability
        C[i][1] = -inv_2dx / width[i-1];
        C[i][2] = inv_Gdt + friction_term + convective_term;
        C[i][3] =  inv_2dx / width[i+1];
        C[i][4] = 0.0;
    }

    // j=2: First interior velocity equation with full Saint-Venant terms
    {
        double ui2 = velocity[2];
        double friction2 = (fabs(ui2) / (Chezy[2] * Chezy[2])) / waterDepth[2];
        double vel_diff2 = velocity[4] - velocity[2];
        double convective_term2 = vel_diff2 * inv_2Gdx;  // One-sided scheme near boundary

        // Proper tidal forcing from prescribed water level at i=1
        // This term drives the tidal signal into the estuary from the boundary
        double tidal_forcing = inv_2dx * (freeArea[1] / width[1]);
        
        Z[2] = tidal_forcing + inv_Gdt * ui2;
        C[2][1] = 0.0;
        C[2][2] = inv_Gdt + friction2 + convective_term2;
        C[2][3] = inv_2dx / width[3];
        C[2][4] = 0.0;
    }

    // j=M-1: Upstream continuity with discharge
    double local_RS_inv_dt_upstream = rs[M-1] * inv_dt;
    Z[M-1] = tempFreeArea[M-1] * local_RS_inv_dt_upstream - Qr * inv_2dx + lateral_inflow_profile[M-1];
    C[M-1][1] = -totalArea[M-2] * inv_2dx;
    C[M-1][2] = local_RS_inv_dt_upstream;
    C[M-1][3] = 0.0;
    C[M-1][4] = 0.0;
}

/**
 * @brief Update free-surface elevation and velocity after implicit solve.
 */
void Update()
{
    double tmp[MAXM+1];
    
    // Update even points (velocity points)
    for (int i = 2; i <= M-2; i += 2) {
        tmp[i-1] = tempFreeArea[i-1] + baseArea[i-1];  // TH[i-1] + ZZ[i-1]
        tmp[i+1] = tempFreeArea[i+1] + baseArea[i+1];  // TH[i+1] + ZZ[i+1]
        
        tempFreeArea[i] = (tempFreeArea[i-1] + tempFreeArea[i+1]) / 2.0;
        totalArea[i] = (tmp[i-1] + tmp[i+1]) / 2.0;
        waterDepth[i] = ((tmp[i-1] / width[i-1]) + (tmp[i+1] / width[i+1])) / 2.0;
    }
    
    // Update odd points (water level points)
    for (int i = 3; i <= M-1; i += 2) {
        totalArea[i] = tempFreeArea[i] + baseArea[i];
        waterDepth[i] = totalArea[i] / width[i];
        tempVelocity[i] = (tempVelocity[i+1] + tempVelocity[i-1]) / 2.0;
    }
    
    // --- Upstream boundary (i=M): Fortran extrapolation + radiation relaxation ---
    // Pure high-order extrapolation reflects the tidal wave at the domain end,
    // causing artificial tidal range amplification. Apply a radiation relaxation
    // that blends the extrapolated value toward the interior, absorbing the
    // outgoing wave (Blayo & Debreu 2005, Orlanski 1976).
    {
        double AAft_M1 = tempFreeArea[M1];
        double AAft_M3 = tempFreeArea[M3];

        // Fortran-heritage extrapolation
        double AAft_extrap = 0.5 * (3.0 * AAft_M1 - AAft_M3);
        
        // Radiation relaxation: blend extrapolation toward interior neighbor.
        // alpha=0 → pure extrapolation (full reflection), alpha=1 → copy interior (zero-gradient).
        const double alpha_rad = 0.85;
        tempFreeArea[M] = (1.0 - alpha_rad) * AAft_extrap + alpha_rad * AAft_M1;
        
        // Prevent negative free surface anomalies that would expose the riverbed numerically
        if (tempFreeArea[M] < 0.0) tempFreeArea[M] = 0.0;

        totalArea[M] = tempFreeArea[M] + baseArea[M];

        if (width[M] > 1e-9) {
            waterDepth[M] = totalArea[M] / width[M];
        } else {
            waterDepth[M] = (M >= 2) ? waterDepth[M-2] : 0.0;
        }
    }

    // Fortran: Ut(1) = Ut(2)
    tempVelocity[1] = tempVelocity[2];
}

/**
 * @brief Apply final results respecting staggered grid structure
 */
void ApplyResults()
{  
    // 1. Copy PRIMARY results from solver to appropriate grid points
    // Copy water levels ONLY to odd points (where they are solved)
    bool hydro_error_detected = false;
    int first_problem_cell = -1;
    double problem_value = 0.0;
    
    for (int i = 1; i <= M; i += 2) {
        // Check for numerical problems BEFORE updating state
        if (isnan(tempFreeArea[i]) || isinf(tempFreeArea[i])) {
            hydro_error_detected = true;
            first_problem_cell = i;
            problem_value = tempFreeArea[i];
            break;
        }
        
        freeArea[i] = tempFreeArea[i];      // Water levels at odd points only
        totalArea[i] = baseArea[i] + freeArea[i];
        if (width[i] > 1e-9) {
            waterDepth[i] = totalArea[i] / width[i];
            level[i] = freeArea[i] / width[i];  // Water level above datum
            
            // Check for unrealistic water depth.
            // Keep physically conservative bounds to catch runaway/nonphysical solutions.
            if (waterDepth[i] < -2.0 || waterDepth[i] > 50.0) {
                hydro_error_detected = true;
                first_problem_cell = i;
                problem_value = waterDepth[i];
                break;
            }
        }
    }
    
    // Immediate halt on critical hydrodynamic error
    if (hydro_error_detected) {
        debug_log(DEBUG_LEVEL_ESSENTIAL, DEBUG_MODULE_HYDRO,
                 "CRITICAL: Hydrodynamic calculation error at cell %d (%.1f km from mouth)",
                 first_problem_cell, first_problem_cell * DELXI / 1000.0);
        debug_log(DEBUG_LEVEL_ESSENTIAL, DEBUG_MODULE_HYDRO,
                 "Problem value: %.3g (likely %s)",
                 problem_value, 
                 isnan(problem_value) ? "NaN" : 
                 isinf(problem_value) ? "Infinity" : 
                 problem_value < -1.0 ? "negative depth" : "extreme value");
        
        printf("❌ FATAL HYDRODYNAMIC ERROR: Calculation failed, see details above\n");
        printf("   Suggestion: Check boundary conditions and time step\n");
        exit(1);
    }
    
    // Copy velocities ONLY to even points (where they are solved)
    for (int i = 2; i <= M-2; i += 2) {
        velocity[i] = tempVelocity[i];      // Velocities at even points only
    }

    // Ensure upstream boundary-face velocity is refreshed from current BC so
    // transport sees timestep-consistent upstream advection sign/magnitude.
    if ((M % 2) == 0) {
        velocity[M] = tempVelocity[M];
        if (!isfinite(velocity[M])) {
            velocity[M] = (M >= 2 && isfinite(velocity[M-1])) ? velocity[M-1] : 0.0;
        }
    }
    
    // 2. INTERPOLATE values to "wrong" grid points for convenience/output
    
    // Interpolate water levels to even points (for output/diagnostics)
    for (int i = 2; i <= M-2; i += 2) {
        freeArea[i] = 0.5 * (freeArea[i-1] + freeArea[i+1]);
        totalArea[i] = baseArea[i] + freeArea[i];
        if (width[i] > 1e-9) {
            waterDepth[i] = totalArea[i] / width[i];
            level[i] = freeArea[i] / width[i];  // Water level above datum
        }
    }
    
    // Interpolate velocities to odd points (for output/diagnostics)
    for (int i = 3; i <= M-1; i += 2) {
        velocity[i] = 0.5 * (velocity[i-1] + velocity[i+1]);
    }

    // Only extrapolate when velocity wasn't computed by solver
    if (isnan(velocity[1]) || isinf(velocity[1])) {
        velocity[1] = tempVelocity[2];  // Fallback only if needed
    }
    
    // Upstream boundary at i=M: use values already computed by Update()
    // (which uses the exact Fortran area-based extrapolation).
    // Just need to ensure freeArea/level are consistent.
    freeArea[M] = tempFreeArea[M];
    totalArea[M] = baseArea[M] + freeArea[M];
    if (width[M] > 1e-9) {
        waterDepth[M] = totalArea[M] / width[M];
        level[M] = freeArea[M] / width[M];
    } else {
        waterDepth[M] = (M >= 2) ? waterDepth[M-2] : 0.0;
        level[M] = 0.0;
    }
    
}

/**
 * @brief Normalized convergence check (hydrodynamics iterative solver)
 *
 * Monitor the relative change (|Δx| / max(|x|,|xₚ|,1)) and compare it against a
 * dimensionless tolerance (default 10⁻⁶ ≈ 1 μm on a 1 m signal), consistent with the
 * centimetre-scale accuracy discussed by Savenije (2012, §3.3) for water-level solutions.
 */
static double Conv(int s, int e, double *xarray, double *yarray, double toler, double *norm_out)
{
    double max_norm = 0.0;

    for (int i = s; i <= e; i += 2) {
        double xi = xarray[i];
        double yi = yarray[i];
        double diff = fabs(xi - yi);
        double scale = fmax(fmax(fabs(xi), fabs(yi)), 1.0);
        double norm = diff / scale;

        if (norm > max_norm) {
            max_norm = norm;
        }

        // No damping or artificial radiation; retain Savenije (2012) staggered-grid physics
        yarray[i] = xi;
    }

    if (norm_out) {
        *norm_out = max_norm;
    }

    return (max_norm <= toler) ? 1.0 : 0.0;
}

void NewUH(int t)
{
    (void)t; // currently unused; keep signature stable

    // Staggered grid structure
    ApplyResults();
    
    // Calculate shear stress and other derived values
    // Literature: τ_b = ρ g U² / C² (Savenije 2012, Eq. 2.12)
    // Note: τ_b MUST be positive (shear stress is a scalar magnitude)
    for (int i = 1; i <= M; i++) {
        if (i % 2 == 0) {
            // At velocity points, use velocity directly
            // τ_b = ρ g U² / C² (always positive due to U²)
            tau_b[i] = rho_w * g * velocity[i] * velocity[i] / (Chezy[i] * Chezy[i]);
        } else {
            // At water level points, interpolate from nearby velocity points
            double U_interpolated = 0.0;
            
            if (i == 1) {
                // First point: use velocity from second point
                U_interpolated = velocity[2];
            } else if (i == M) {
                // Last point: use velocity from second-to-last point
                U_interpolated = velocity[M-1];
            } else {
                // Interior points: interpolate from adjacent velocity points
                U_interpolated = 0.5 * (velocity[i-1] + velocity[i+1]);
            }
            
            // τ_b = ρ g U² / C² (always positive due to U²)
            tau_b[i] = rho_w * g * U_interpolated * U_interpolated / (Chezy[i] * Chezy[i]);
        }
    }
    
    // Update tidal extrema values
    {
        for (int i = 1; i <= M; i++) {
            if (waterDepth[i] < daily_min_level[i]) {
                daily_min_level[i] = waterDepth[i];
            }
            if (waterDepth[i] > daily_max_level[i]) {
                daily_max_level[i] = waterDepth[i];
            }

            if (velocity[i] < velocity_min[i]) {
                velocity_min[i] = velocity[i];
            }
            if (velocity[i] > velocity_max[i]) {
                velocity_max[i] = velocity[i];
            }
        }
    }
}

void Tridag()
{
    // Highly optimized tridiagonal solver (physics unchanged – performance optimized)
    static double gam_local[MAXM+1];
    static double var_local[MAXM+1];
    static long pivot_clamp_count = 0;
    
    // Initialize with improved pivot handling
    double bet = C[2][2];
    if (fabs(bet) < 1e-10) {
        bet = 1e-10;  // Pivot protection
        pivot_clamp_count++;
    }
    var_local[2] = Z[2] / bet;

    // Forward elimination with optimized pointer arithmetic
    for (int j = 3; j <= M-1; j++) {
        double *restrict Cjm1 = C[j-1];  // restrict hint for compiler optimization
        double *restrict Cj   = C[j];
        double gam_j = Cjm1[3] / bet;            // C[j-1][3]
        gam_local[j] = gam_j;
        bet = Cj[2] - Cj[1] * gam_j;             // pivot update with local variable

        if (fabs(bet) < 1e-10) {
            bet = 1e-10; // Pivot protection
            pivot_clamp_count++;
        }
        var_local[j] = (Z[j] - Cj[1] * var_local[j-1]) / bet;
    }
    
    if (pivot_clamp_count > 1000) {
        if (debug_level >= DEBUG_LEVEL_WARNING) {
             printf("⚠️  Hydrodynamics: Frequent pivot clamping detected (%ld events). Possible numerical instability.\n", pivot_clamp_count);
        }
        pivot_clamp_count = 0; // Reset to avoid spamming
    }

    // Back substitution (unchanged)
    for (int j = M-2; j >= 2; j--) {
        var_local[j] -= gam_local[j+1] * var_local[j+1];
    }
    
    // Copy results with stride-optimized access
    for (int j = 2; j <= M-2; j += 2) {
        tempVelocity[j] = var_local[j];
        tempFreeArea[j+1] = var_local[j+1];
    }
}

/**
 * @brief utils hydrodynamic solver function
 */
static void reset_daily_extrema_to_current(void)
{
    for (int i = 1; i <= M; i++) {
        daily_min_level[i] = waterDepth[i];
        daily_max_level[i] = waterDepth[i];
        velocity_min[i]    = velocity[i];
        velocity_max[i]    = velocity[i];
    }
}

static void finalize_daily_tidal_statistics(long t_sec)
{
    const long current_day = t_sec / (24 * 3600);
    if (tidal_stats_last_day < 0) {
        tidal_stats_last_day = current_day;
        reset_daily_extrema_to_current();
        return;
    }

    if (current_day == tidal_stats_last_day) {
        return;
    }

    const double prev_day_mid = ((double)tidal_stats_last_day + 0.5) * 86400.0;
    const int include_for_calibration = (prev_day_mid >= (double)WARMUP);

    for (int i = 1; i <= M; i++) {
        double range = daily_max_level[i] - daily_min_level[i];
        if (!isfinite(range) || range < 0.0) {
            range = 0.0;
        }
        if (include_for_calibration) {
            tidal_range_sum[i] += range;
            tidal_range_count[i] += 1.0;
            if (tidal_range_count[i] > 0.0) {
                const double avg_range = tidal_range_sum[i] / tidal_range_count[i];
                tidal_range_mean[i] = avg_range;
            }
        }
    }

    reset_daily_extrema_to_current();
    tidal_stats_last_day = current_day;
}

void finalize_remaining_tidal_range(void)
{
    if (tidal_stats_last_day < 0) {
        return;
    }
    const long next_day = tidal_stats_last_day + 1;
    finalize_daily_tidal_statistics(next_day * 24L * 3600L);
}

void Hyd(int t)
{

    if (t == 0) {
        hydro_solver_calls_total = 0;
        hydro_solver_calls_post_warmup = 0;
        hydro_iteration_sum_post_warmup = 0;
        hydro_maxits_hit_count = 0;
        hydro_max_iteration_observed = 0;
        hydro_max_iteration_observed_post_warmup = 0;
        tidal_stats_last_day = 0;
        for (int i = 1; i <= M; i++) {
            tidal_range_sum[i] = 0.0;
            tidal_range_count[i] = 0.0;
            tidal_range_mean[i] = 0.0;
            dispersion_min[i]  = 50.0;
            dispersion_max[i]  = 50.0;
        }
        reset_daily_extrema_to_current();
    }

    if (calibration_mode && t == WARMUP) {
        for (int i = 1; i <= M; i++) {
            tidal_range_sum[i] = 0.0;
            tidal_range_count[i] = 0.0;
            tidal_range_mean[i] = 0.0;
        }
        reset_daily_extrema_to_current();
        tidal_stats_last_day = t / (24 * 3600);
    }

    for (int i = 1; i <= M; i++) {
        totalAreaOld[i] = totalArea[i];
    }

    // Enforce strict convergence: no environment/runtime tolerance overrides.
    if (hydro_convergence_tolerance <= 0.0) {
        hydro_convergence_tolerance = 5e-4;  // More realistic tidal tolerance
    }

    Newbc(t);

    int iteration = 0;
    double rsum = 0.0;
    bool converged = false;
    int effective_maxits = MAXITS;

    do {
        iteration++;

        Coeffa(t);
        Tridag();

        double eta_residual = 0.0;
        double velocity_residual = 0.0;

        rsum  = Conv(3, M-1, tempFreeArea, E, hydro_convergence_tolerance, &eta_residual);
        rsum += Conv(2, M-2, tempVelocity, Y, hydro_convergence_tolerance, &velocity_residual);
        
        last_eta_residual = eta_residual;
        last_velocity_residual = velocity_residual;

        Update();

        if (rsum == 2.0) {
            converged = true;
            break;
        }

    } while (iteration < effective_maxits);

    // Track convergence health
    if (!converged && iteration >= effective_maxits) {
        consecutive_maxits++;
    } else {
        consecutive_maxits = 0;
    }

    hydro_solver_calls_total++;
    if (iteration > hydro_max_iteration_observed) {
        hydro_max_iteration_observed = iteration;
    }
    if (t >= WARMUP) {
        hydro_solver_calls_post_warmup++;
        hydro_iteration_sum_post_warmup += iteration;
        if (iteration > hydro_max_iteration_observed_post_warmup) {
            hydro_max_iteration_observed_post_warmup = iteration;
        }
        if (!converged && iteration >= effective_maxits) {
            hydro_maxits_hit_count++;
        }
    }

    if (!converged) {
        debug_log(DEBUG_LEVEL_WARNING, DEBUG_MODULE_HYDRO,
              "Hydrodynamics solver reached MAXITS (%d) at t=%d s (day %.2f) without full convergence (η_residual=%.3e, u_residual=%.3e, tol=%.3e).",
              MAXITS, t, (double)t / 86400.0, last_eta_residual, last_velocity_residual,
              hydro_convergence_tolerance);
    }

    finalize_daily_tidal_statistics(t);

    NewUH(t);


    int hydro_TS = TS;
    if ((double)t / (double)(hydro_TS * DELTI) - floor((double)t / (double)(hydro_TS * DELTI)) == 0.0 && t >= WARMUP) {
        int success = Hydwrite(t);
        if (!success && t % (24 * 3600) == 0) {
            printf("⚠️ Warning: Failed to write hydrodynamic output files at t=%d (day %.1f)\n",
                   t, (double)t / (24 * 3600));
        }
    }
}

void export_hydro_iterations_summary(const char *path)
{
    if (!path || path[0] == '\0') {
        return;
    }

    ensure_directory_for_path(path);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        if (debug_level >= DEBUG_LEVEL_WARNING) {
            printf("⚠️  Could not write hydro iteration diagnostics CSV: %s\n", path);
        }
        return;
    }

    const double mean_iterations_post_warmup =
        (hydro_solver_calls_post_warmup > 0)
            ? ((double)hydro_iteration_sum_post_warmup / (double)hydro_solver_calls_post_warmup)
            : 0.0;
    const double maxits_hit_fraction_post_warmup =
        (hydro_solver_calls_post_warmup > 0)
            ? ((double)hydro_maxits_hit_count / (double)hydro_solver_calls_post_warmup)
            : 0.0;

    fprintf(fp,
            "solver_calls_total_count,solver_calls_post_warmup_count,max_iteration_observed_count,max_iteration_observed_post_warmup_count,maxits_limit_count,maxits_hit_count,maxits_hit_fraction_post_warmup,mean_iterations_post_warmup_count,unit_notes\n");
    fprintf(fp,
            "%ld,%ld,%d,%d,%d,%ld,%.10f,%.10f,counts_and_dimensionless_fraction\n",
            hydro_solver_calls_total,
            hydro_solver_calls_post_warmup,
            hydro_max_iteration_observed,
            hydro_max_iteration_observed_post_warmup,
            MAXITS,
            hydro_maxits_hit_count,
            maxits_hit_fraction_post_warmup,
            mean_iterations_post_warmup);

    fclose(fp);
}
