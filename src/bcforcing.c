/**
 * @file bcforcing.c
 * @brief Implementation of boundary conditions for hydrodynamic calculations
 *
 */

#include "define.h" // Include global definitions
#include "variables.h" // Include global variable declarations
#include "diagnostics.h" // Consolidated diagnostics system
#include <stddef.h>  // For NULL
#include <math.h>
#include <stdbool.h>
#include <stdarg.h>  // For variable argument functions (va_start, va_end)
#include <string.h>

// Forward declarations for configuration-driven smoothing controls (params.txt)
extern int enable_upstream_discharge_smoothing;
extern double upstream_discharge_smoothing_timescale_days;

// The hardcoded is_urban_tributary_name() has been removed.
// We strictly use the configuration-driven Tributary->is_urban flag now.

// Modified Bessel function I0 approximation (Numerical Recipes/Cephes style)
// Max relative error ~1e-7, sufficient for forcing normalization.
static double bessel_i0_approx(double x) {
    double ax = fabs(x);
    if (ax < 3.75) {
        double y = x / 3.75;
        y *= y;
        return 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 +
               y * (0.2659732 + y * (0.0360768 + y * 0.0045813)))));
    }
    double y = 3.75 / ax;
    return (exp(ax) / sqrt(ax)) *
           (0.39894228 + y * (0.01328592 + y * (0.00225319 +
           y * (-0.00157565 + y * (0.00916281 + y * (-0.02057706 +
           y * (0.02635537 + y * (-0.01647633 + y * 0.00392377))))))));
}

static double apply_upstream_discharge_smoothing(int t, double raw_discharge) {
    static double smoothed_value = 0.0;
    static int has_state = 0;
    static int last_update_time = -1;

    if (!enable_upstream_discharge_smoothing || upstream_discharge_smoothing_timescale_days <= 0.0) {
        smoothed_value = raw_discharge;
        has_state = 1;
        last_update_time = t;
        return raw_discharge;
    }

    double tau_seconds = upstream_discharge_smoothing_timescale_days * 86400.0;
    if (!(tau_seconds > 0.0) || isnan(tau_seconds)) {
        smoothed_value = raw_discharge;
        has_state = 1;
        last_update_time = t;
        return raw_discharge;
    }

    // Reset state when simulation restarts or rewinds
    if (!has_state || t < last_update_time) {
        smoothed_value = raw_discharge;
        has_state = 1;
        last_update_time = t;
        return smoothed_value;
    }

    if (t == last_update_time) {
        return smoothed_value;
    }

    double dt = (double)(t - last_update_time);
    if (dt < 0.0) {
        dt = 0.0;
    }

    double alpha = 1.0 - exp(-dt / tau_seconds);
    if (alpha < 0.0) {
        alpha = 0.0;
    } else if (alpha > 1.0) {
        alpha = 1.0;
    }

    smoothed_value += alpha * (raw_discharge - smoothed_value);
    last_update_time = t;
    return smoothed_value;
}

/**
 * @brief Unified function to interpolate input time series data with better logging.
 * 
 * @param t Current simulation time (seconds)
 * @param timeArray Array of time values (days since start)
 * @param dataArray Array of corresponding data values
 * @param dataSize Size of the data arrays
 * @param isHourly Whether data is hourly (true) or daily (false)
 * @param dataName Name of the data for error messages
 * @param defaultValue Default value to use if data is missing
 * @return Interpolated value at time t
 */
double interpolateInputData(int t, double* timeArray, double* dataArray, int dataSize, 
                           bool isHourly, const char* dataName, double defaultValue)
{    
    // Validate pointers and size
    if (!dataArray || dataSize <= 0) {
        printf("❌ ERROR: Invalid data array for %s (size=%d)\n", dataName, dataSize);
        return defaultValue;
    }

    // Hourly series (e.g., tidal elevation, light) MUST remain dynamic during warmup.
    // For analysis, outputs begin at t=WARMUP; to avoid a forcing phase discontinuity,
    // align hourly forcings to the same (t-WARMUP) clock used by daily series.
    if (isHourly) {
        // Use a phase-shifted, wrapped clock during warmup; use (t-WARMUP) thereafter.
        const double warmup_hours = (double)WARMUP / 3600.0;
        double hours_since_start;
        if (t < WARMUP) {
            // Shift warmup so that it ends exactly at hour 0 of the analysis record.
            hours_since_start = ((double)t / 3600.0) - warmup_hours;
        } else {
            hours_since_start = (double)(t - WARMUP) / 3600.0;
        }

        // If timeArray exists, wrap by its total span and interpolate *periodically*
        // (including interpolation across last->first) to ensure continuity at wrap points.
        if (timeArray) {
            const double series_start = timeArray[0];
            double spacing = (dataSize > 1) ? fabs(timeArray[1] - timeArray[0]) : 0.0;
            if (spacing <= 0.0) {
                spacing = 1.0;
            }
            const double total_span = fmax(spacing, (timeArray[dataSize - 1] - series_start) + spacing);

            double wrapped = fmod(hours_since_start, total_span);
            if (wrapped < 0.0) wrapped += total_span;
            const double target = series_start + wrapped;

            // Binary search bracket: time[lo] <= target < time[hi]
            int lo = 0, hi = dataSize - 1;
            if (target <= timeArray[0]) return dataArray[0];

            // Periodic interpolation across the wrap boundary.
            if (target >= timeArray[dataSize - 1]) {
                const double t_lo = timeArray[dataSize - 1];
                const double t_hi = series_start + total_span;
                const double dt = t_hi - t_lo;
                const double f = (dt > 0.0) ? (target - t_lo) / dt : 0.0;
                const double val = dataArray[dataSize - 1] + f * (dataArray[0] - dataArray[dataSize - 1]);
                return (isnan(val) || isinf(val)) ? defaultValue : val;
            }

            while (hi - lo > 1) {
                int mid = (lo + hi) >> 1;
                if (timeArray[mid] <= target) lo = mid; else hi = mid;
            }
            double dt = timeArray[hi] - timeArray[lo];
            double f = (dt > 0.0) ? (target - timeArray[lo]) / dt : 0.0;
            f = fmax(0.0, fmin(1.0, f));
            double val = dataArray[lo] + f * (dataArray[hi] - dataArray[lo]);
            return (isnan(val) || isinf(val)) ? defaultValue : val;
        } else {
            // Unit-spaced hours; wrap by length
            double period_end = (double)(dataSize - 1);
            if (period_end < 1.0) return dataArray[0];
            double target = fmod(hours_since_start, period_end);
            if (target < 0.0) target += period_end;
            int idx = (int)floor(target);
            if (idx < 0) return dataArray[0];
            if (idx >= dataSize - 1) return dataArray[dataSize - 1];
            double f = target - idx;
            double val = dataArray[idx] + f * (dataArray[idx + 1] - dataArray[idx]);
            return (isnan(val) || isinf(val)) ? defaultValue : val;
        }
    }

    if (dataSize <= 0) {
        return defaultValue;
    }

    // When the daily series does not carry an explicit time axis, treat it as
    // uniformly spaced and wrap the warmup period over the available samples.
    if (!timeArray) {
        if (dataSize == 1) {
            return (isnan(dataArray[0]) || isinf(dataArray[0])) ? defaultValue : dataArray[0];
        }

        if (t < WARMUP) {
            double warmup_days = (double)t / 86400.0;
            double total_span = (double)dataSize;
            if (total_span <= 0.0) total_span = 1.0;
            double warmup_length_days = (double)WARMUP / 86400.0;
            double phase = warmup_days - warmup_length_days;
            double wrapped = fmod(phase, total_span);
            if (wrapped < 0.0) wrapped += total_span;
            int idx = (int)floor(wrapped);
            double f = wrapped - idx;
            int next_idx = (idx + 1) % dataSize;
            double val = dataArray[idx] + f * (dataArray[next_idx] - dataArray[idx]);
            return (isnan(val) || isinf(val)) ? defaultValue : val;
        }

        double target_days = (double)(t - WARMUP) / 86400.0;
        int idx = (int)floor(target_days);
        if (idx < 0) return dataArray[0];
        if (idx >= dataSize - 1) return dataArray[dataSize - 1];
        double f = target_days - idx;
        double val = dataArray[idx] + f * (dataArray[idx + 1] - dataArray[idx]);
        return (isnan(val) || isinf(val)) ? defaultValue : val;
    }

    // Daily series with explicit time axis: treat the series as periodic over its sampled span.
    // During warmup we wrap to reach dynamic equilibrium; after warmup we keep the same
    // wrapped clock to avoid dataset-end discontinuities and to match the forcing-looping
    // behavior used historically in biogeochemical helpers.
    double series_start = timeArray[0];
    double spacing = (dataSize > 1) ? fabs(timeArray[1] - timeArray[0]) : 0.0;
    if (spacing <= 0.0) {
        spacing = 1.0; // Fall back to daily cadence when metadata absent or non-monotonic
    }
    double total_span = fmax(spacing, (timeArray[dataSize - 1] - series_start) + spacing);

    double target_days;
    if (t < WARMUP) {
        const double warmup_days = (double)t / 86400.0;
        const double warmup_length_days = (double)WARMUP / 86400.0;
        const double phase = warmup_days - warmup_length_days;
        double wrapped = fmod(phase, total_span);
        if (wrapped < 0.0) wrapped += total_span;
        target_days = series_start + wrapped;
    } else {
        const double sim_days = (double)(t - WARMUP) / 86400.0;
        double wrapped = fmod(sim_days, total_span);
        if (wrapped < 0.0) wrapped += total_span;
        target_days = series_start + wrapped;
    }

    // Binary search for daily series using provided time array
    int lo = 0, hi = dataSize - 1;
    if (target_days <= timeArray[0]) return dataArray[0];
    // Periodic interpolation across the wrap boundary.
    if (target_days >= timeArray[dataSize - 1]) {
        const double t_lo = timeArray[dataSize - 1];
        const double t_hi = series_start + total_span;
        const double dt = t_hi - t_lo;
        const double f = (dt > 0.0) ? (target_days - t_lo) / dt : 0.0;
        const double val = dataArray[dataSize - 1] + f * (dataArray[0] - dataArray[dataSize - 1]);
        return (isnan(val) || isinf(val)) ? defaultValue : val;
    }
    while (hi - lo > 1) {
        int mid = (lo + hi) >> 1;
        if (timeArray[mid] <= target_days) lo = mid; else hi = mid;
    }
    double dt = timeArray[hi] - timeArray[lo];
    double f = (dt > 0.0) ? (target_days - timeArray[lo]) / dt : 0.0;
    f = fmax(0.0, fmin(1.0, f));
    double val = dataArray[lo] + f * (dataArray[hi] - dataArray[lo]);
    return (isnan(val) || isinf(val)) ? defaultValue : val;
}

/**
 * @brief Computes the tidal elevation at time @p t using physics-based multi-component tidal forcing.
 * 
 * This implementation is based on Savenije's theory of tidal propagation in estuaries with:
 * - M2 (Principal lunar semidiurnal): Main tidal component (period ~12.42 hours)
 * - S2 (Principal solar semidiurnal): Secondary tidal component (period 12.00 hours)
 * - M4 (Shallow water overtide of M2): Important for tidal asymmetry (period ~6.21 hours)
 */
double Tide(int t)
{    
    // Add caching to prevent redundant calculations
    static int last_tide_t = -1;
    static double cached_tide = 0.0;
    
    // If we've already calculated tide for this timestep, return cached value
    if (t == last_tide_t) {
        return cached_tide;
    }
    
    // Update cache timestamp
    last_tide_t = t;
    
    // Use actual boundary elevation data if available
    static int elevation_stats_initialized = 0;
    static double elevation_mean = 0.0;

    if (forcingData[FORCING_ELEVATION].dataSize > 0) {
        if (!elevation_stats_initialized) {
            int n = forcingData[FORCING_ELEVATION].dataSize;
            if (n > 0) {
                double sum = 0.0;
                for (int i = 0; i < n; ++i) {
                    sum += forcingData[FORCING_ELEVATION].data[i];
                }
                elevation_mean = sum / (double)n;
            } else {
                elevation_mean = 0.0;
            }
            elevation_stats_initialized = 1;
        }
        // Interpolate from the hourly elevation data
        double eta = interpolateInputData(t,
                                        forcingData[FORCING_ELEVATION].time,
                                        forcingData[FORCING_ELEVATION].data,
                                        forcingData[FORCING_ELEVATION].dataSize,
                                        true, // Allow extrapolation for safety
                                        "Elevation",
                                        0.0); // Default value if interpolation fails
        
        // Apply amplification factor around the dataset mean (Savenije 2012 tidal attenuation calibration)
        double eta_centered = eta - elevation_mean;
        eta = elevation_mean + eta_centered * AMPL_CORRECTION;

        // Store in cache and return scaled CSV data
        cached_tide = eta;
        return cached_tide;
    }
    
    // Physics-based synthetic tide following Savenije's theory
    // ⚠️ WARNING: This fallback mode activates when no elevation CSV data is available
    static int synthetic_mode_warned = 0;
    if (!synthetic_mode_warned) {
        printf("⚠️  WARNING: No elevation CSV data found. Using SYNTHETIC TIDE mode.\n");
        printf("   → For production runs, provide hourly elevation data in config_input.txt\n");
        synthetic_mode_warned = 1;
    }
    
    // Define the mean sea level and tidal parameters
    static double MEAN_SEA_LEVEL = 0.0;      // Will be set based on AMPL
    
    static double AMPL_M2, AMPL_S2;          // Main tidal components
    static double PERIOD_M2 = 0.0;           // M2 tidal period - from params.txt
    static double PERIOD_S2 = 0.0;           // S2 tidal period - from params.txt
    
    // Phase relationships following Savenije's research on tidal asymmetry
    static double PHI_M2, PHI_S2, PHI_M4;    // Phases for physical realism
    static double omega_M2, omega_S2;        // Angular frequencies
    
    
    static int initialized = 0;
    
    // Initialize amplitudes and phases based on AMPL parameter and physical understanding
    // following Savenije's research on tidal dynamics in estuaries
    if (initialized == 0) {
        // Load parameters from configuration (Directive 2: Configuration-driven)
        PERIOD_M2 = TIDAL_PERIOD_M2_HOURS;     // From params.txt
        PERIOD_S2 = TIDAL_PERIOD_S2_HOURS;     // From params.txt
        
        MEAN_SEA_LEVEL = 0.0; // Mean sea level at datum - CORRECTED
        
        AMPL_M2 = AMPL * M2_fraction;  // Configurable M2 component fraction (AMPL is amplitude)
        AMPL_S2 = AMPL * S2_fraction;  // Configurable S2 component fraction (AMPL is amplitude)
        
        // Initialize phases with physically consistent relationships from configuration
        // The phase relationships determine tidal asymmetry characteristics
        PHI_M2 = TIDAL_PHASE_M2_RAD;     // From params.txt (Directive 2)
        PHI_S2 = TIDAL_PHASE_S2_RAD;     // From params.txt (Directive 2)
        
        // M4 phase relationship determines whether there's flood or ebb dominance
        PHI_M4 = 2.0 * PHI_M2;         // Standard M4-M2 relationship
        
        // Angular frequencies directly from tidal theory
        omega_M2 = 2.0 * M_PI / (PERIOD_M2 * 3600.0);
        omega_S2 = 2.0 * M_PI / (PERIOD_S2 * 3600.0);
        
        initialized = 1;
        
        printf("🌊 Physics-based tide model initialized following Savenije (2012):\n");
        printf("   • AMPL=%.2f m, MSL=%.2f m, M2=%.2f m, S2=%.2f m\n", 
               AMPL, MEAN_SEA_LEVEL, AMPL_M2, AMPL_S2);
        printf("   • Tidal period: %.2f hours (M2), %.2f hours (S2)\n",
               PERIOD_M2, PERIOD_S2);
        printf("   • Tidal range at mouth: %.2f m\n", 2.0 * AMPL);
        printf("   • M4/M2 amplitude ratio: %.2f (creates tidal asymmetry)\n", M4_ratio);
        printf("   • Phase configuration: %s dominant\n", 
               (PHI_M4 < 0.5) ? "flood" : "ebb");
    }
    
    // Calculate primary astronomical constituents (M2, S2)
    double eta = MEAN_SEA_LEVEL + 
                AMPL_M2 * sin(omega_M2 * t + PHI_M2) +
                AMPL_S2 * sin(omega_S2 * t + PHI_S2);
    
    // Add M4 shallow water overtide (essential for tidal asymmetry)
    // The M4 component is generated by nonlinear interaction of M2 with itself
    double A_M4 = AMPL_M2 * M4_ratio;    // M4 amplitude proportional to M2
    double phi_M4 = 2.0 * PHI_M2 + PHI_M4; // Phase relationship following tidal theory
    eta += A_M4 * sin(2.0 * omega_M2 * t + phi_M4);

    // Apply configuration-driven amplification (Savenije 2012 tidal damping calibration)
    double tidal_component = eta - MEAN_SEA_LEVEL;
    eta = MEAN_SEA_LEVEL + tidal_component * AMPL_CORRECTION;
    
    // 🌊 STEADY-STATE WARMUP: Keep FULL tidal dynamics during warmup
    
    // Store in cache and return
    cached_tide = eta;
    return cached_tide;
}

/**
 * @brief  Returns upstream NET discharge for time @p t.
 * @param  t  Simulation time in seconds.
 * @return     Net upstream discharge (m³/s).
 *
 * @details
 *  - Negative sign indicates flow direction into the model domain (per grid convention).
 *  - Linear interpolation is performed between known daily discharge data points.
 *  - Optional exponential smoothing can be applied via enable_upstream_discharge_smoothing.
 *  - During warmup: Uses FULL physics (no artificial ramping) to ensure proper equilibration.
 *  - Validation checks prevent unrealistic discharge values from destabilizing the model.
 *
 * Scientific basis: Savenije (2012) §2.3 recommends smoothing coarse discharge records
 * to avoid spurious high-frequency oscillations in the salinity response.
 */
double Discharge_ups(int t) //Freshwater inflow (upstream), This is the net discharge, NOT instant discharge!
{
    static double prev_discharge = -100.0; // Track previous discharge for smoothing
    double default_discharge = -100.0;     // Safe default value, negative for inflow
    
    // Calculate the mean discharge from the first few days of data to use as reference
    static double mean_discharge = 0.0;
    static int mean_initialized = 0;
    
    if (mean_initialized == 0 && forcingData[FORCING_DISCHARGE].dataSize > 0 && forcingData[FORCING_DISCHARGE].data) {
        // Calculate mean from first few days or all data if less available
        int sample_size = fmin(10, forcingData[FORCING_DISCHARGE].dataSize);
        double sum = 0.0;
        for (int i = 0; i < sample_size; i++) {
            sum += forcingData[FORCING_DISCHARGE].data[i];
        }
        mean_discharge = sum / sample_size;
        mean_initialized = 1;
        default_discharge = mean_discharge;
    }

    // Guard: missing discharge forcing data.
    // Avoid calling interpolateInputData() with empty arrays (which prints a hard error).
    if (!forcingData[FORCING_DISCHARGE].data || !forcingData[FORCING_DISCHARGE].time || forcingData[FORCING_DISCHARGE].dataSize <= 0) {
        static int warned_missing = 0;
        if (!warned_missing && debug_level >= DEBUG_LEVEL_WARNING) {
            printf("⚠️ Upstream discharge forcing missing/empty; using default discharge %.3g m³/s\n", default_discharge);
            warned_missing = 1;
        }
        return default_discharge;
    }
    
    // Get discharge from measured data with negative sign for inflow
    // Same treatment during warmup and post-warmup: FULL physics, no artificial ramping
    // This ensures proper hydrodynamic equilibration (Blumberg & Mellor 1987)
    double discharge = -interpolateInputData(t,
                               forcingData[FORCING_DISCHARGE].time,
                               forcingData[FORCING_DISCHARGE].data,
                               forcingData[FORCING_DISCHARGE].dataSize,
                               false,  // Discharge data is daily
                               "upstream discharge",
                               default_discharge);

    // Optional sensitivity scaling (dimensionless). Applied before smoothing so
    // the filtered series corresponds to the scaled forcing.
    if (isfinite(upstream_discharge_scale) && upstream_discharge_scale > 0.0 && upstream_discharge_scale != 1.0) {
        discharge *= upstream_discharge_scale;
    }
    
    // Apply optional low-pass filtering (Savenije 2012 §2.3)
    discharge = apply_upstream_discharge_smoothing(t, discharge);

    prev_discharge = discharge;
    return discharge;
}

/**
 * @brief Mass-conserving physics-based tributary discharge calculation with tidal modulation
 * @param t Current simulation time in seconds
 * @param tributary_index Index of the tributary in the tributaries array
 * @return Modulated tributary discharge [m³/s] (positive = inflow to estuary)
 * 
 * @details
 * This function implements a MASS-CONSERVING tidal modulation of tributary discharge
 * that accounts for:
 * 1. Backwater effects from tidal stage variations
 * 2. Tidal pumping mechanisms (flood vs ebb tide effects)
 * 3. Mass conservation over tidal cycles (daily average preserved)
 */
/**
 * @brief Get tributary discharge from CSV data - clean scientific approach
 * @details Simply interpolates the CSV discharge data without arbitrary modifications.
 */
double calculatePhysicsBasedTributaryDischarge(int t, int tributary_index) {
    if (tributary_index < 0 || tributary_index >= numTributaries) {
        return 0.0;
    }
    
    // Validate tributary information
    int cell = tributaries[tributary_index].cellIndex;
    if (cell < 1 || cell > M || cell >= MAXM) return 0.0;
    
    // Get discharge from CSV data - no modifications, just interpolation
    int dataSize = tributaries[tributary_index].dischargeDataSize;
    if (dataSize <= 0) return 0.0;
    
    double discharge = interpolateInputData(
                            t, 
                            tributaries[tributary_index].dischargeTime,
                            tributaries[tributary_index].discharge,
                            dataSize,
                            false,  // daily data
                            "tributary discharge",
                            0.0);
    
    if (!(discharge > 0.0)) {
        return 0.0;
    }

    // Optional subdaily pulse modulation (mean-preserving over a full cycle):
    //   q(t) = Q_daily * exp(a*sin(theta)) / I0(a)
    // where mean_theta[exp(a*sin(theta))] = I0(a), so E[q]=Q_daily.
    if (enable_tributary_subdaily_pulses && tributary_pulse_shape > 0.0) {
        double a = tributary_pulse_shape;
        if (tributaries[tributary_index].is_urban) {
            a *= urban_tributary_pulse_boost;
        }
        if (!isfinite(a) || a < 0.0) {
            a = 0.0;
        }
        if (a > 4.0) {
            a = 4.0; // prevent unrealistically sharp spikes
        }

        if (a > 0.0) {
            const double omega_M2 = 2.0 * M_PI / (TIDAL_PERIOD_M2_HOURS * 3600.0);
            const double theta = omega_M2 * (double)t;
            const double norm = bessel_i0_approx(a);
            if (norm > 1e-12 && isfinite(norm)) {
                const double mod = exp(a * sin(theta)) / norm;
                if (isfinite(mod) && mod > 0.0) {
                    discharge *= mod;
                }
            }
        }
    }

    return (discharge > 0.0 && isfinite(discharge)) ? discharge : 0.0;
}

/**
 * @brief  Returns discharge at any given point and time, with option to get tributary contribution.
 * @param  t  Simulation time (seconds).
 * @param  i  Grid index (1 to M, where 1 is downstream mouth, M is upstream).
 * @param  trib_discharge  Optional pointer to store tributary discharge component (pass NULL if not needed).
 *                         When provided, the absolute value of tributary discharge (positive) will be stored.
 * @return     Discharge at the specified grid point (m³/s).
 *
 * @details
 *  This function computes the discharge at any point in the model domain:
 *  - For upstream boundary (i=M): calls Discharge_ups()
 *  - For downstream boundary (i=1): accumulates all tributary flows + upstream flow
 *  - For interior points: accumulates only tributary flows that enter at or upstream of point i
 *  - When tributaries are disabled, discharge remains constant from upstream to downstream
 */
double Discharge(int t, int i, double *trib_discharge)
{
    // Fischer (1979) mode: tributaries handled as mass sources in transport equation
    // Return only the base upstream discharge - no tributary modification
    if (tributary_mass_injection_mode == TRIB_INJECTION_FISCHER) {
        double Q_up = Discharge_ups(t);
        if (trib_discharge) {
            *trib_discharge = 0.0; // Tributaries handled in transport, not discharge
        }
        return Q_up;
    }
    
    // Original C implementation for RUTHERFORD mode
    // Return NET freshwater discharge magnitude (≥ 0) at cross-section i
    // Upstream NET discharge is negative (into domain) in Discharge_ups(t)
    // Tributary discharges are always positive
    static int last_t = -1;
    static double cached_Q_up_mag = 0.0;

    if (t != last_t) {
        cached_Q_up_mag = fabs(Discharge_ups(t)); // magnitude only
        last_t = t;
    }

    // Bounds and trivial cases
    if (i >= M) {
        if (trib_discharge) *trib_discharge = 0.0;
        return cached_Q_up_mag; // upstream section: only main river NET freshwater
    }
    if (i < 1) {
        if (trib_discharge) *trib_discharge = 0.0;
        return cached_Q_up_mag; // clamp
    }

    double Q_trib_accumulated = 0.0;

    if (tributaryEnabled && tributaries && numTributaries > 0) {
        for (int j = 0; j < numTributaries; j++) {
            // Accumulate tributaries located at or downstream of cross-section i (toward mouth)
            if (tributaries[j].cellIndex >= i) {
                double qj = calculatePhysicsBasedTributaryDischarge(t, j);
                if (qj > 0.0 && !isnan(qj) && !isinf(qj)) {
                    Q_trib_accumulated += qj; // always positive
                }
            }
        }
    }

    if (trib_discharge) *trib_discharge = Q_trib_accumulated;

    // NET freshwater magnitude (non-negative) along space
    return cached_Q_up_mag + Q_trib_accumulated;
}


/**
 * @brief Updates biogeochemical boundary conditions with flow-direction awareness
 * @param t Current simulation time (s)
 * 
 * This function updates boundary condition values (clb and cub) but respects
 * transport-calculated mixing for salinity to allow dispersion effects to work.
 */
void bgboundary(int t)
{
    int var;
    double prev_clb, prev_cub;
    
    /* Process each biogeochemical variable */
    for (var = 0; var < CHEM_COUNT; var++) {
        /* Skip auto-calculated variables (don't need boundary conditions) */
        if (var == pCO2 || var == PH || var == CO2) continue;
        
        /* Store previous values for detecting large jumps */
        prev_clb = v[var].clb;
        prev_cub = v[var].cub;
        
        if (t < WARMUP) {
            // WARMUP: Use first data value for steady state equilibration
            if (upstreamBC[var].time && upstreamBC[var].data && upstreamBC[var].dataSize > 0) {
                v[var].cub = upstreamBC[var].data[0];
            }
            
            if (downstreamBC[var].time && downstreamBC[var].data && downstreamBC[var].dataSize > 0) {
                v[var].clb = downstreamBC[var].data[0];
            }
        }
        else {
            // MAIN SIMULATION: Apply CSV boundary data directly (configuration-driven)
            // Following Memory Bank Directive 2: All parameters from configuration files
            // CRITICAL: Prevent model state contamination of boundary conditions
            
            // Set upstream boundary condition from CSV data with strict validation
            if (upstreamBC[var].time && upstreamBC[var].data && upstreamBC[var].dataSize > 0) {
                double new_cub = interpolateInputData(
                    t,
                    upstreamBC[var].time,
                    upstreamBC[var].data,
                    upstreamBC[var].dataSize,
                    false,  // daily boundary data
                    variableNames[var],
                    prev_cub
                );
                
                // Enhanced validation to prevent boundary condition corruption
                // Only accept reasonable boundary values to prevent mass balance violations
                bool is_valid_boundary = true;
                if (isnan(new_cub) || isinf(new_cub)) {
                    is_valid_boundary = false;
                } else if (new_cub < 0.0) {
                    // Negative concentrations violate physical constraints
                    is_valid_boundary = false;
                } else if (var == Sal && (new_cub > 40.0 || new_cub < 0.0)) {
                    // Salinity must be physically reasonable (0-40 ppt)
                    is_valid_boundary = false;
                } else if (var != Sal && new_cub > 1e6) {
                    // Other variables: prevent extremely large values that could destabilize model
                    is_valid_boundary = false;
                }
                
                if (!is_valid_boundary) {
                    fprintf(stderr, "\n❌ FATAL ERROR: Invalid upstream boundary data for %s at day %.2f (value=%.3f)\n", 
                           variableNames[var], (double)(t - WARMUP)/86400.0, new_cub);
                    fprintf(stderr, "❌ Simulation aborted to prevent non-mechanistic results. Please check input CSV files.\n");
                    exit(EXIT_FAILURE);
                }
                
                v[var].cub = new_cub;
            }
            
            // Set downstream boundary condition from CSV data with strict validation
            if (downstreamBC[var].time && downstreamBC[var].data && downstreamBC[var].dataSize > 0) {
                double new_clb = interpolateInputData(
                    t,
                    downstreamBC[var].time,
                    downstreamBC[var].data,
                    downstreamBC[var].dataSize,
                    false,  // daily boundary data
                    variableNames[var],
                    prev_clb
                );
                
                // Enhanced validation to prevent boundary condition corruption
                bool is_valid_boundary = true;
                if (isnan(new_clb) || isinf(new_clb)) {
                    is_valid_boundary = false;
                } else if (new_clb < 0.0) {
                    // Negative concentrations violate physical constraints
                    is_valid_boundary = false;
                } else if (var == Sal && (new_clb > 40.0 || new_clb < 0.0)) {
                    // Salinity must be physically reasonable (0-40 ppt)
                    is_valid_boundary = false;
                } else if (var != Sal && new_clb > 1e6) {
                    // Other variables: prevent extremely large values that could destabilize model
                    is_valid_boundary = false;
                }
                
                if (!is_valid_boundary) {
                    fprintf(stderr, "\n❌ FATAL ERROR: Invalid downstream boundary data for %s at day %.2f (value=%.3f)\n", 
                           variableNames[var], (double)(t - WARMUP)/86400.0, new_clb);
                    fprintf(stderr, "❌ Simulation aborted to prevent non-mechanistic results. Please check input CSV files.\n");
                    exit(EXIT_FAILURE);
                }
                
                v[var].clb = new_clb;
        
            }
        }
                
        // Reset flux arrays to ensure clean state
        for (int i = 0; i <= M; i++) {
            v[var].concflux[i] = 0.0;
            v[var].advflux[i] = 0.0;
            v[var].disflux[i] = 0.0;
        }
    }

    /*------------------------------------------------------*/
    /* ✅ Process Tributaries - Configuration-driven approach */
    /*------------------------------------------------------*/
    if (tributaryEnabled) {
        for (int j = 0; j < numTributaries; j++) {
            // Update tributary concentrations from CSV data for each variable
            for (int var = 0; var < CHEM_COUNT; var++) {
                // Skip auto-calculated variables
                if (var == pCO2 || var == PH || var == CO2) continue;
                
                // Get CSV data arrays
                double *timeArray = tributaries[j].chemicalData[var].timeArray;
                double *dataArray = tributaries[j].chemicalData[var].dataArray;
                int dataSize = tributaries[j].chemicalData[var].dataSize;
                
                if (timeArray && dataArray && dataSize > 0) {
                    if (t < WARMUP) {
                        // During warmup: use first value for steady state
                        tributaries[j].concentration[var] = dataArray[0];
                    } else {
                        // Main simulation: interpolate CSV data directly
                        double trib_conc = interpolateInputData(t,
                                                              timeArray,
                                                              dataArray,
                                                              dataSize,
                                                              false,  // daily data
                                                              tributaries[j].name,
                                                              tributaries[j].concentration[var]);
                        
                        // Only check for data corruption (NaN/Inf/negative), not "jumps"
                        // Natural concentration variations in CSV data should be preserved
                        if (isnan(trib_conc) || isinf(trib_conc) || trib_conc < 0.0) {
                            fprintf(stderr, "\n❌ FATAL ERROR: Invalid tributary concentration for %s-%s at day %.2f (value=%.3f)\n", 
                                   tributaries[j].name, variableNames[var], (double)(t-WARMUP)/86400.0, trib_conc);
                            fprintf(stderr, "❌ Simulation aborted to prevent non-mechanistic results. Please check input CSV files.\n");
                            exit(EXIT_FAILURE);
                        } else {
                            tributaries[j].concentration[var] = trib_conc;
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Calculates wind speed at time t by interpolating from daily wind data.
 * 
 * @param t Current simulation time in seconds
 * @return Wind speed value interpolated for the current time
 */
double WS(int t)
{
    double ws_value, restday, day;
    int nday;

    day = floor((double)t / (24.0 * 60.0 * 60.0));
    restday = ((double)t / (24.0 * 60.0 * 60.0)) - day;
    nday = (int)day;

    // Check if we have wind data from forcingData structure
    if (!forcingData[FORCING_WIND].data) {
        return 0.0; // Safety check
    }

    if (day < 0.0)
    {
        ws_value = forcingData[FORCING_WIND].data[0];
    }
    else
    {
        // Bounds check using actual data size
        if (nday >= 0 && nday+1 < forcingData[FORCING_WIND].dataSize) {
            ws_value = (forcingData[FORCING_WIND].data[nday] + 
                      (forcingData[FORCING_WIND].data[nday + 1] - forcingData[FORCING_WIND].data[nday]) * restday);
        } else if (nday >= 0 && nday < forcingData[FORCING_WIND].dataSize) {
            ws_value = forcingData[FORCING_WIND].data[nday]; // Last valid entry
        } else {
            ws_value = forcingData[FORCING_WIND].data[0]; // Default to first value
        }
    }

    return ws_value;
}

/**
 * @brief  Computes wind speed at time @p t using an exponential decay upstream.
 * @param  t  Simulation time (s).
 * @param  i  Grid index along the channel.
 * @return    Adjusted wind speed (m/s) at the cell @p i.
 */
double windspeed(int t, int i)
{
    double default_wind = 0.0;
    
    // Get base wind speed
    double Wspeed = interpolateInputData(t,
                                       forcingData[FORCING_WIND].time,
                                       forcingData[FORCING_WIND].data,
                                       forcingData[FORCING_WIND].dataSize,
                                       false,  // Wind data is daily
                                       "wind speed",
                                       default_wind);
    
    // Apply spatial decay
    return Wspeed * exp(-(((i + 1) - 1) * ((double)DELXI)) / ((double)EL));
}

/**
 * @brief Calculates the tidal phase at time t.
 *
 * @details
 * This function calculates the tidal phase (0-2π) within the tidal cycle.
 * Used to determine whether it's flood or ebb tide, critical for boundary conditions.
 * 
 * @param t Simulation time in seconds.
 * @return Tidal phase in radians (0-2π)
 */
double calculateTidalPhase(int t) {
    // Calculate phase using M2 tidal frequency (dominant component)
    static double omega_M2 = 0.0; // Will be initialized from configuration parameters
    
    // Initialize omega_M2 from configuration on first call (Directive 2: Configuration-driven)
    if (omega_M2 == 0.0) {
        omega_M2 = 2.0 * M_PI / (TIDAL_PERIOD_M2_HOURS * 3600.0); // From params.txt
    }
    
    // Extract phase from sine wave (returns 0 to 2π)
    double phase = fmod(omega_M2 * t, 2.0 * M_PI);
    if (phase < 0) phase += 2.0 * M_PI;  // Ensure positive phase

    return phase;
}