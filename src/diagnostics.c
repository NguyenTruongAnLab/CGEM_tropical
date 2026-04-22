/**
 * @file diagnostics.c
 * @brief Diagnostics, validation, and debugging system for C-GEM model
 * @author Nguyen Truong An
 * @date Last updated: 11/2025
 * 
 * Provides:
 * - Debug logging with configurable levels
 * - Output control for calibration mode suppression
 * - Numerical stability checks
 * - Mass computation utilities
 */

#include "diagnostics.h"
#include "define.h"
#include "variables.h"
#include "file.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <float.h>

// External variables
extern double disp[MAXM + 1];
extern double totalArea[MAXM + 1];
extern double velocity[MAXM + 1];
extern int M;
extern int DELTI;
extern int DELXI;
extern long WARMUP;
extern double level[MAXM + 1];
extern int debug_level;
extern int calibration_mode;

// Global output control state (merged from output_control.c)
int global_output_enabled = 1;

// Debug system globals
static time_t last_debug_output_time[10] = {0};
static int last_reported_milestone_day = -1;

// Module names for debug output
static const char* module_names[] = {
    "GENERAL", "INIT", "HYDRO", "TRANSPORT", "BC", "BIOGEO", "FILE", "COUPLING"
};

// Coupling diagnostics
static HydroTransportCouplingData coupling_data = {0, 0.0, 0.0, 0.0, false, false, false, 0};

// Validation summary (minimal)
static ValidationSummary validation_summary = {0};

// Calibration-safe fatal latch: when calibration_mode>0, do not terminate the process.
// Instead, latch a fatal request so the calibration driver can penalize the parameter set.
static bool g_fatal_requested = false;

// =====================================================================
// DIAGNOSTICS INITIALIZATION
// =====================================================================

void diagnostics_init(int level) {
    debug_level = level;
    memset(&validation_summary, 0, sizeof(validation_summary));
    validation_summary.min_salt_intrusion_km = 1e9;
    validation_summary.min_dispersion_m2s = 1e9;
    init_coupling_diagnostics();
}

// =====================================================================
// DEBUG LOGGING SYSTEM
// =====================================================================

bool debug_is_enabled(DebugLevel level) {
    return (debug_level >= level);
}

const char* debug_get_module_name(DebugModule module) {
    if (module < DEBUG_MODULE_COUNT) {
        return module_names[module];
    }
    return "UNKNOWN";
}

static bool should_debug_output(DebugLevel level, DebugModule module) {
    if (!debug_is_enabled(level)) return false;
    if (calibration_mode > 0) return false;  // Silent during calibration
    if (level == DEBUG_LEVEL_ESSENTIAL) return true;
    
    int current_day = (int)((time(NULL) - WARMUP) / 86400.0);
    static int last_debug_day[DEBUG_MODULE_COUNT] = {0};
    
    if (level == DEBUG_LEVEL_WARNING) {
        time_t now = time(NULL);
        if (last_debug_output_time[module] == 0 || current_day > last_debug_day[module]) {
            last_debug_output_time[module] = now;
            last_debug_day[module] = current_day;
            return true;
        }
        return false;
    }
    
    if (level == DEBUG_LEVEL_DETAIL) {
        if (is_milestone_day(time(NULL)) || last_debug_output_time[module] == 0) {
            time_t now = time(NULL);
            if (now - last_debug_output_time[module] > DEBUG_FREQ_DAILY) {
                last_debug_output_time[module] = now;
                return true;
            }
        }
        return false;
    }
    return false;
}

void debug_log(DebugLevel level, DebugModule module, const char* format, ...) {
    if (!should_debug_output(level, module)) return;
    
    va_list args;
    va_start(args, format);
    
    const char* severity = (level == DEBUG_LEVEL_WARNING) ? "⚠️ WARNING: " :
                           (level == DEBUG_LEVEL_DETAIL) ? "🔍 " : "ℹ️ INFO: ";
    
    printf("[%s] %s", debug_get_module_name(module), severity);
    vprintf(format, args);
    printf("\n");
    fflush(stdout);
    
    va_end(args);
}

void debug_log_periodic(DebugLevel level, DebugModule module, 
                        int frequency_seconds, const char *key, 
                        const char *format, ...) {
    if (!debug_is_enabled(level)) return;
    
    static time_t last_periodic_times[100] = {0};
    static char periodic_keys[100][64];
    static int num_periodic_keys = 0;
    
    int key_index = -1;
    for (int i = 0; i < num_periodic_keys; i++) {
        if (strcmp(periodic_keys[i], key) == 0) {
            key_index = i;
            break;
        }
    }
    
    if (key_index == -1 && num_periodic_keys < 100) {
        key_index = num_periodic_keys;
        strncpy(periodic_keys[key_index], key, 63);
        periodic_keys[key_index][63] = '\0';
        num_periodic_keys++;
    }
    
    if (key_index == -1) return;
    
    time_t now = time(NULL);
    if (now - last_periodic_times[key_index] < frequency_seconds) return;
    last_periodic_times[key_index] = now;
    
    va_list args;
    va_start(args, format);
    printf("[%s] ", debug_get_module_name(module));
    vprintf(format, args);
    printf("\n");
    fflush(stdout);
    va_end(args);
}

bool is_milestone_day(int t) {
    int current_day = (int)((t - WARMUP) / 86400.0);
    
    if (current_day != last_reported_milestone_day) {
        if (current_day == 1 || current_day == 30 || current_day == 90 || 
            current_day == 180 || current_day == 365 || current_day % 365 == 0) {
            last_reported_milestone_day = current_day;
            return true;
        }
    }
    return false;
}

// =====================================================================
// OUTPUT CONTROL (merged from output_control.c)
// =====================================================================

void init_output_control(int enabled) {

    global_output_enabled = enabled;
}

int printf_controlled(const char *format, ...) {
    if (calibration_mode > 0 || !global_output_enabled) return 0;
    va_list args;
    va_start(args, format);
    int result = vprintf(format, args);
    va_end(args);
    return result;
}

void set_output_enabled(int enabled) { global_output_enabled = enabled; }
int is_output_enabled(void) { return global_output_enabled; }

void init_coupling_diagnostics(void) {
    memset(&coupling_data, 0, sizeof(coupling_data));
    coupling_data.last_update_time = -1;
}

// =====================================================================
// NUMERICAL STABILITY CHECK
// =====================================================================

static void fatal_tracer_error(const char* stage, int t, int s, int i, double value, const char* reason)
{
    // In calibration mode, latch the first fatal only and allow the caller to unwind.
    if (calibration_mode > 0) {
        if (g_fatal_requested) {
            return;
        }
        g_fatal_requested = true;
    }

    // Calibration performance: avoid flooding stderr in optimization runs.
    // Users can raise debug_level (>=WARNING) to see the detailed fatal messages.
    if (calibration_mode > 0 && debug_level < DEBUG_LEVEL_WARNING) {
        return;
    }

    const char* name = getSpeciesName(s);
    if (!name) name = "UNKNOWN";
    // Stable ASCII prefix for test harness regex matching.
    fprintf(stderr,
            "FATAL TRACER ERROR: stage=%s t=%d species=%d (%s) cell=%d value=%.17g reason=%s\n",
            stage ? stage : "(null)",
            t,
            s,
            name,
            i,
            value,
            reason ? reason : "(null)");
    fflush(stderr);
    if (calibration_mode > 0) {
        return;
    }
    exit(EXIT_FAILURE);
}

void diagnostics_fatal_tracer_error(const char* stage, int t, int s, int i, double value, const char* reason)
{
    fatal_tracer_error(stage, t, s, i, value, reason);
}

void diagnostics_validate_tracer_state(const char* stage, int t)
{
    if (g_fatal_requested) {
        return;
    }

    // Tolerance for numerical precision near zero (floating point artifacts).
    // In highly eutrophic systems with rapid nutrient depletion (e.g., strong
    // phytoplankton uptake depleting NO3 near zero, or P adsorption onto SPM),
    // the implicit reaction step can produce values slightly below zero.
    // These are numerical artifacts, not physical, and should be clamped to zero.
    const double NEG_TOL = -1e-4;

    // Enforce invariant on all prognostic tracers (exclude purely diagnostic carbonate outputs).
    for (int s = 0; s < CHEM_COUNT; ++s) {
        if (v[s].env != 1) {
            continue;
        }
        if (s == pCO2 || s == PH || s == CO2) {
            continue;
        }
        for (int i = 1; i <= M; ++i) {
            double val = v[s].c[i];
            if (!isfinite(val)) {
                fatal_tracer_error(stage, t, s, i, val, "nonfinite");
                return;
            }
            // Clamp tiny negatives (floating point artifacts) to zero
            if (val < 0.0 && val >= NEG_TOL) {
                v[s].c[i] = 0.0;
                val = 0.0;
            }
            if (val < NEG_TOL) {
                fatal_tracer_error(stage, t, s, i, val, "negative");
                return;
            }
        }
    }
}


bool diagnostics_fatal_requested(void)
{
    return g_fatal_requested;
}

void diagnostics_clear_fatal_request(void)
{
    g_fatal_requested = false;
}

bool check_numerical_stability(int t) {
    static int last_check_time = 0;
    static int problem_count = 0;
    
    if (t - last_check_time < 60) return true;
    last_check_time = t;
    // Physical sanity thresholds (generous bounds for any estuary)
    static const double MAX_VELOCITY    = 15.0;    // [m/s]
    static const double MIN_DEPTH       = 0.01;    // [m]
    static const double MAX_DEPTH       = 100.0;   // [m]
    static const double MAX_DISPERSION  = 50000.0; // [m²/s]

    bool has_nan = false;
    bool has_extreme = false;
    int problem_loc = -1;
    char problem_field[32] = "";
    double extreme_val = 0.0;
    
    for (int i = 1; i <= M; i++) {
        if (!isfinite(velocity[i])) {
            has_nan = true; problem_loc = i;
            strcpy(problem_field, "velocity"); break;
        }
        if (fabs(velocity[i]) > MAX_VELOCITY) {
            has_extreme = true; problem_loc = i;
            extreme_val = velocity[i];
            strcpy(problem_field, "velocity"); break;
        }
        
        if (!isfinite(waterDepth[i])) {
            has_nan = true; problem_loc = i;
            strcpy(problem_field, "waterDepth"); break;
        }
        if (waterDepth[i] < MIN_DEPTH || waterDepth[i] > MAX_DEPTH) {
            has_extreme = true; problem_loc = i;
            extreme_val = waterDepth[i];
            strcpy(problem_field, "waterDepth"); break;
        }
        
        if (!isfinite(disp[i])) {
            has_nan = true; problem_loc = i;
            strcpy(problem_field, "dispersion"); break;
        }
        if (disp[i] > MAX_DISPERSION) {
            has_extreme = true; problem_loc = i;
            extreme_val = disp[i];
            strcpy(problem_field, "dispersion"); break;
        }
    }
    
    if (has_nan) {
        debug_log(DEBUG_LEVEL_ESSENTIAL, DEBUG_MODULE_GENERAL,
                 "CRITICAL: %s[%d] contains NaN/Inf at t=%d", problem_field, problem_loc, t);
        return false;
    }
    
    if (has_extreme) {
        problem_count++;
        if (problem_count <= 3) {
            debug_log(DEBUG_LEVEL_WARNING, DEBUG_MODULE_GENERAL,
                     "Extreme value: %s[%d] = %.3g at t=%d", problem_field, problem_loc, extreme_val, t);
        }
        if (problem_count >= 5) {
            debug_log(DEBUG_LEVEL_ESSENTIAL, DEBUG_MODULE_GENERAL,
                     "CRITICAL: Persistent extreme values - instability detected");
            return false;
        }
    } else {
        if (problem_count > 0) problem_count--;
    }
    
    return true;
}
