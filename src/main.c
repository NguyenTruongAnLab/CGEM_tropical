/**
 * @file main.c
 * @brief Main entry point and control flow for the C-GEM (Carbon-Generic Estuarine Model).
 * @author An Nguyen
 * @date Last updated: 11/2025
 *
 * C-GEM is a 1D reactive-transport model for estuarine systems, simulating
 * hydrodynamics, sediment transport, and biogeochemical processes.
 */

#include "define.h"
#include "variables.h"
#include "utilities.h"
#include "biogeo.h"
#include "diagnostics.h"
#include "transport.h"
#include "init.h"
#include "file.h"

#ifndef CGEM_ENABLE_CALIBRATION
#define CGEM_ENABLE_CALIBRATION 0
#endif

#if CGEM_ENABLE_CALIBRATION
#include "calibration.h"
#endif

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>

int main(int argc, char **argv) {
    (void)argc;
    (void)argv;
    
    // Setup and configuration
    setup_locale_and_signals();

    int parsed_mode = calibration_mode;
    int parsed_debug = debug_level;
    const char *params_path = cgem_get_params_path();
    if (read_calibration_flags(params_path, &parsed_mode, &parsed_debug) != 0) {
        fprintf(stderr, "⚠️ Unable to read calibration flags from %s; using defaults.\n", params_path);
    }
    calibration_mode = parsed_mode;
    debug_level = parsed_debug;

    // Test/CI override: force calibration off regardless of params.txt.
    // This keeps smoke tests deterministic even when researchers keep calibration_mode != 0 locally.
    {
        const char *force_no_cal_env = getenv("CGEM_FORCE_NO_CALIBRATION");
        if (force_no_cal_env && force_no_cal_env[0] != '\0' && atoi(force_no_cal_env) != 0) {
            calibration_mode = 0;
        }
    }

    // Initialize output control system - suppress verbose output during calibration
#if !CGEM_ENABLE_CALIBRATION
    if (calibration_mode != 0) {
        fprintf(stderr,
                "\u26a0\ufe0f Calibration requested (calibration_mode=%d) but CGEM_ENABLE_CALIBRATION=0; disabling calibration and continuing normal simulation.\n",
                calibration_mode);
        calibration_mode = 0;
    }
    init_output_control(1);
#else
    init_output_control(calibration_mode == 0 ? 1 : 0);
#endif

    const char *debug_env = getenv("CGEM_DEBUG_LEVEL");
    if (debug_env && debug_env[0] != '\0') {
        int env_debug = atoi(debug_env);
        if (env_debug >= DEBUG_LEVEL_ESSENTIAL && env_debug <= DEBUG_LEVEL_DETAIL) {
            debug_level = env_debug;
        }
    }

#if CGEM_ENABLE_CALIBRATION
    if (calibration_mode != 0) {
        fprintf(stderr, "[CALIB] Entering calibration (mode=%d stage=%d) params=%s\n",
                calibration_mode, calibration_stage, cgem_get_params_path());
        fflush(stderr);
        int calibration_status = run_hydrodynamic_calibration_from_params();
        fprintf(stderr, "[CALIB] Calibration returned status=%d\n", calibration_status);
        fflush(stderr);
        return (calibration_status == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
    }
#endif

    // Track simulation timing
    clock_t start_time = clock();

    // Initialize model
    Init();

    // Runtime override: transport-only mode (no reactions, no erosion/deposition)
    // This is intended for boundary/lateral validation workflows.
    {
        const char *transport_only_env = getenv("CGEM_TRANSPORT_ONLY");
        if (transport_only_env && transport_only_env[0] != '\0' && atoi(transport_only_env) != 0) {
            enable_biogeochemical_reactions = 0;
            enable_suspended_sediment_dynamics = 0;
            // Keep carbonate diagnostics on by default so pH/pCO2/CO2 are consistent with DIC/AT transport.
            enable_carbonate_diagnostics = 1;
        }

        // Always print the effective mode (helps avoid confusion with validators that only read params.txt)
        if (debug_level >= DEBUG_LEVEL_ESSENTIAL) {
            printf("ℹ️ Run mode flags: enable_biogeochemical_reactions=%d, enable_carbonate_diagnostics=%d, enable_suspended_sediment_dynamics=%d\n",
                   enable_biogeochemical_reactions, enable_carbonate_diagnostics, enable_suspended_sediment_dynamics);
        }
    }

    // Validate and override WARMUP/MAXT via environment variables
    const char *warmup_env = getenv("CGEM_WARMUP_DAYS"); // Warmup timeframe 
    if (warmup_env && warmup_env[0] != '\0') {
        long warmup_days = (long)atof(warmup_env);
        if (warmup_days >= 0) {
            WARMUP = warmup_days * 86400;
        }
    }

    const char *maxt_env = getenv("CGEM_MAXT_DAYS");    // Total sim time (warmup + sim)
    if (maxt_env && maxt_env[0] != '\0') {
        long maxt_days = (long)atof(maxt_env);
        if (maxt_days > 0) {
            MAXT = maxt_days * 86400;
        }
    }

    // Optional diagnostics override: shorten main simulation duration via environment variable (scientific debugging only)
    const char *mass_diag_days_env = getenv("CGEM_MASS_DIAG_DAYS");
    if (mass_diag_days_env && mass_diag_days_env[0] != '\0') {
        double diag_days = atof(mass_diag_days_env);
        if (diag_days > 0.0) {
            long warmup_steps = WARMUP / DELTI;
            long diagnostic_steps = (long)((diag_days * 86400.0) / (double)DELTI);
            long new_total_steps = warmup_steps + diagnostic_steps;
            long new_MAXT = new_total_steps * DELTI;
            if (new_MAXT < MAXT) {
                if (debug_level >= DEBUG_LEVEL_DETAIL) {
                    printf("🔬 CGEM_MASS_DIAG_DAYS=%g ➜ Reducing MAXT from %ld to %ld seconds (warmup + %.1f days)\n",
                           diag_days, MAXT, new_MAXT, diag_days);
                }
                MAXT = new_MAXT;
            }
        }
    }

    if (debug_level >= DEBUG_LEVEL_ESSENTIAL) {
        printf("ℹ️ Effective runtime controls: warmup=%.2f d, total=%.2f d, dt=%d s\n",
               (double)WARMUP / 86400.0,
               (double)MAXT / 86400.0,
               DELTI);
    }
    
    // Initialize simulation variables
    long steps_done = 0;
    double last_printed_day = -5.0;
    double progress_pct = 0;
    clock_t last_progress_time = start_time;
    long t_start = 0;
    
    initialize_simulation_variables(start_time, &steps_done, &last_printed_day,
                                   &progress_pct, &last_progress_time, &t_start);

    // =====================================================================
    // MAIN SIMULATION LOOP
    // =====================================================================
    
    for (long t = t_start; t <= MAXT; t += DELTI) {
        steps_done++;
        Hyd(t);
        bgboundary(t);
        Transport(t);

        if (t > 0 && !check_numerical_stability(t)) {
            printf("\n❌ SIMULATION HALTED: Numerical stability issue at Day %.1f\n", (double)t / (24.0 * 60.0 * 60.0));
            perform_final_cleanup();
            return EXIT_FAILURE;
        }

        if (enable_biogeochemical_reactions || enable_carbonate_diagnostics) {
            Biogeo(t);
        }

        if (enable_suspended_sediment_dynamics) {
            updateSuspendedSediment(t);
        }

        // Write state outputs at the configured interval AFTER all updates.
        // This ensures exported pH/CO2/pCO2 (diagnostics) and eutrophication states
        // reflect the same timestep (no one-step lag).
        {
            const long output_stride = (long)TS * (long)DELTI;
            if (t >= WARMUP && output_stride > 0 && (t % output_stride) == 0) {
                SaveAllWaterQuality((int)t);
                for (int s = 0; s < MAXV; s++) {
                    if (v[s].env == 1 && s != pCO2 && s != PH && s != CO2) {
                        Fluxwrite(s, (int)t);
                    }
                }
            }
        }
    } // End of main simulation loop

    // Export runtime diagnostics required by manuscript reviewer package
    export_hydro_iterations_summary("outputs/metrics/hydro_iterations_summary.csv");
    export_limiter_events_summary("outputs/metrics/limiter_events_summary.csv");
    
    // Final cleanup
    perform_final_cleanup();

    // Post-processing: export to NetCDF (optional).
    // Python validation scripts are separate tools — run them manually after the
    // simulation completes (see README.md and tools/README.md).
    const char *skip_export_env = getenv("CGEM_SKIP_VALIDATION");
    const int skip_export = (skip_export_env && skip_export_env[0] != '\0' && atoi(skip_export_env) != 0);

    if (!skip_export) {
        run_postprocessing_export();
    } else {
        printf("ℹ️ CGEM_SKIP_VALIDATION=1 → Skipping NetCDF export\n");
    }

    const char *output_label = (output_storage_format == OUTPUT_FORMAT_BINARY) ? "BINARY (.bin)" : "CSV";
    printf("✅ All simulation data available as %s files for downstream analysis\n", output_label);

    // Report runtime (useful for I/O benchmarking; physics unaffected)
    {
        double elapsed_sec = (double)(clock() - start_time) / (double)CLOCKS_PER_SEC;
        printf("⏱️  Runtime: %.2f seconds\n", elapsed_sec);
    }

    return 0;
}
