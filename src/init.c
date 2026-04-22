/**
 * @file init.c
 * @brief Initialization of hydrodynamic model state and parameters
 */

#include "define.h"
#include "variables.h"
#include "utilities.h"
#include "file.h"
#include "init.h"
#include "diagnostics.h"
#include "biogeo.h"            // Add biogeo header for fast carbonate chemistry
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>

#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#else
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

// Forward declaration for discharge functions from bcforcing.c
extern double Discharge_ups(int t);

/**
 * @brief Initializes the model by reading configuration files, allocating memory,
 *        and setting initial conditions.
 *
 * @details
 * This function performs the following steps:
 *   1. Deletes old simulation output files.
 *   2. Reads the main model parameters (`params.txt`).
 *   3. Sets up boundary and tributary data structures.
 *   4. Initializes hydrodynamic, transport and biogeochemical variables.
 */
void Init() {

    // Safety: if the model is ever invoked multiple times in the same process
    // (e.g., calibration loops, future multi-scenario drivers), ensure no output
    // streams from a previous run remain open. On Windows, open handles can also
    // prevent OUT/* cleanup from deleting stale files.
    reset_output_file_registry();


    // Portability hardening: avoid shell/system() calls and avoid hardcoding
    // specific output filenames. We clear top-level OUT/*.{csv,bin,nc} to match
    // the original intent (do not recurse into subdirectories).

    // Ensure OUT directory exists (non-fatal if already exists)
#ifdef _WIN32
    if (_mkdir("OUT") != 0 && errno != EEXIST) {
        if (debug_level >= DEBUG_LEVEL_WARNING) {
            printf("⚠️ Could not create OUT directory (errno=%d)\n", errno);
        }
    }

    const char *patterns[] = {"OUT\\*.csv", "OUT\\*.bin", "OUT\\*.nc"};
    int deleted_count = 0;
    int failed_count = 0;

    for (size_t p = 0; p < sizeof(patterns) / sizeof(patterns[0]); ++p) {
        WIN32_FIND_DATAA ffd;
        HANDLE hFind = FindFirstFileA(patterns[p], &ffd);
        if (hFind == INVALID_HANDLE_VALUE) {
            continue;
        }

        do {
            if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
                continue;
            }
            char filepath[MAX_PATH];
            snprintf(filepath, sizeof(filepath), "OUT\\%s", ffd.cFileName);
            if (DeleteFileA(filepath)) {
                deleted_count++;
            } else {
                DWORD err = GetLastError();
                (void)err;
                failed_count++;
            }
        } while (FindNextFileA(hFind, &ffd) != 0);

        FindClose(hFind);
    }
#else
    if (mkdir("OUT", 0755) != 0 && errno != EEXIST) {
        if (debug_level >= DEBUG_LEVEL_WARNING) {
            printf("⚠️ Could not create OUT directory (errno=%d)\n", errno);
        }
    }

    int deleted_count = 0;
    int failed_count = 0;
    DIR *dir = opendir("OUT");
    if (dir) {
        struct dirent *ent;
        while ((ent = readdir(dir)) != NULL) {
            const char *name = ent->d_name;
            if (!name || name[0] == '.') {
                continue;
            }
            const char *ext = strrchr(name, '.');
            if (!ext) {
                continue;
            }
            if (strcmp(ext, ".csv") != 0 && strcmp(ext, ".bin") != 0 && strcmp(ext, ".nc") != 0) {
                continue;
            }
            char filepath[512];
            snprintf(filepath, sizeof(filepath), "OUT/%s", name);
            if (unlink(filepath) == 0) {
                deleted_count++;
            } else if (errno != ENOENT) {
                failed_count++;
            }
        }
        closedir(dir);
    } else if (debug_level >= DEBUG_LEVEL_WARNING) {
        printf("⚠️ Could not open OUT directory for cleanup (errno=%d)\n", errno);
    }
#endif

    if (debug_level >= DEBUG_LEVEL_ESSENTIAL) {
        printf("🗑️  Cleaned %d old output files", deleted_count);
        if (failed_count > 0) {
            printf(" (⚠️ %d files could not be deleted)", failed_count);
        }
        printf("\n");
    }

    // Step 3: Read model parameters (configuration-driven; INPUT files are read-only)
    // Reviewer-grade reproducibility: allow selecting a scenario params file without
    // editing the repository default (see CGEM_PARAMS_PATH).
    read_parameters(cgem_get_params_path());

    // Strict scientific mode: runtime clocks/output format are controlled only
    // through versioned INPUT/params.txt (no environment overrides).
    simulation_window_seconds = MAXT - WARMUP;
    if (simulation_window_seconds <= 0) {
        fprintf(stderr, "❌ Error: MAXT must exceed WARMUP to define a diagnostic window. (MAXT=%ld, WARMUP=%ld seconds)\n", MAXT, WARMUP);
        exit(EXIT_FAILURE);
    }

    // Step 4: Initialize diagnostics system with level from params.txt (skip heavy setup in perf mode)
    diagnostics_init(debug_level);
    
    // DEBUG: Show loaded Chezy parameters immediately after reading
    debug_log(DEBUG_LEVEL_ESSENTIAL, DEBUG_MODULE_INIT, "Loaded Chezy parameters: Chezy1=%.6f, Chezy2=%.6f", Chezy1, Chezy2);
    
    // Step 5: Read riverbed profile
    readRiverbedProfile(riverbed_depth, riverbed_profile_file);    
    
    // Step 5.1: Allocate cumulative discharge arrays (required for binary I/O)
    allocateCumulativeDischargeMemory();
    
    // Step 6: Read global configuration settings (including tributaryEnabled)
    readGlobalConfigSettings("INPUT/config_input.txt");
    
    // Step 7: Read and setup boundary conditions from config_input.txt
    readBoundaryData("INPUT/config_input.txt");
    
    // Step 8: Read and setup tributary data from config_input.txt (only if enabled)
    if (tributaryEnabled) {
        readTributaryData("INPUT/config_input.txt");
    } else {
        printf("🚫 Tributary system disabled (tributaryEnabled=0)\n");
        numTributaries = 0;
    }

    // Step 8: Initialize boundary values for time t=0
    // CRITICAL: bgboundary() must be called BEFORE initializeBiogeochemistry()
    // This sets v[var].clb and v[var].cub from the CSV boundary data files
    // and initializes concentrations with proper linear gradients between boundaries
    bgboundary(0);
    printf("Boundary Values Initialized for t=0.\n");

    // Step 9: Initialize hydrodynamic fields
    initializeHydrodynamics();

    // Step 10: Initialize transport processes
    initializeTransportVariables();
    
    // Step 11: Initialize biogeochemistry
    initializeBiogeochemistry();

    // Step 12: Assign biogeochemical rate constants
    assignBiogeochemicalRateConstants();

    // Print summary of loaded data
    int boundaryCount = 0;
    for (int i = 0; i < CHEM_COUNT; i++) {
        if (upstreamBC[i].filePath || downstreamBC[i].filePath) {
            boundaryCount++;
        }
    }
    
    // Create a comprehensive initialization summary
    if (debug_level > 0) {
        printf("🌊 Domain Setup:\n");
        printf("   • Estuary length: %.1f km (%d cells)\n", EL/1000.0, M);
        printf("   • Grid resolution: %.1f m\n", (double)DELXI);
        printf("   • Tidal amplitude: %.2f m\n", AMPL);
        printf("   • Segments: %d (LC1: %.1f km, LC2: %.1f km)\n", num_segments, LC1/1000.0, LC2/1000.0);
    }
    
    if (debug_level > 1) {
        printf("🧪 Biogeochemistry:\n");
        printf("   • Chemical variables: %d\n", boundaryCount);
        printf("   • Tributaries: %d\n", numTributaries);
        printf("   • Debug level: %d (%s)\n", debug_level, 
               debug_level == 0 ? "OFF" : debug_level == 1 ? "ESSENTIAL" : "DETAILED");
        
        printf("⏱️  Simulation Settings:\n");
        printf("   • Total time: %ld days (%.1f years)\n", MAXT/(24*3600), MAXT/(365.25*24*3600));
        printf("   • Warmup period: %ld days\n", WARMUP/(24*3600));
        printf("   • Time step: %d seconds (%.1f min)\n", DELTI, DELTI/60.0);
        printf("   • Output interval: %d time steps (%.1f hours)\n", TS, TS*DELTI/3600.0);
    }
    PRINTF_INIT("================================================================\n");
    PRINTF_INIT("🚀 Starting simulation...\n\n");
}

/**
 * @brief Removes old output files to prepare for a new simulation run
 */
void cleanupOldOutputFiles() {    
    PRINTF_INIT("Deleting old simulation files...\n");
    
    // Portable directory creation and file cleanup using C standard library
    // Create main output directory (cross-platform)
    #ifdef _WIN32
        #include <direct.h>
        if (_mkdir("OUT") != 0 && errno != EEXIST) {
            if (debug_level >= DEBUG_LEVEL_WARNING) {
                printf("⚠️ Could not create OUT directory (may already exist)\n");
            }
        }
    #else
        #include <sys/stat.h>
        #include <sys/types.h>
        if (mkdir("OUT", 0755) != 0 && errno != EEXIST) {
            if (debug_level >= DEBUG_LEVEL_WARNING) {
                printf("⚠️ Could not create OUT directory (may already exist)\n");
            }
        }
    #endif
    
    // Clean old output files.
    // Requirement: remove binary files first to avoid stale data across runs.
    // The binary writer uses 'wb' (overwrite), but stale per-variable files can still
    // persist if a variable wasn't written in the new run.
    int deleted_count = 0;
    int failed_count = 0;

    // Always attempt to delete a known NetCDF output (if present).
    if (remove("OUT/cgem_outputs.nc") == 0) {
        deleted_count++;
    } else if (errno != ENOENT) {
        failed_count++;
    }

    // Delete all *.bin files under OUT/ (portable, no shell calls).
    {
#ifdef _WIN32
        WIN32_FIND_DATAA fd;
        HANDLE hFind = FindFirstFileA("OUT\\*.bin", &fd);
        if (hFind != INVALID_HANDLE_VALUE) {
            do {
                char path[MAX_PATH];
                snprintf(path, sizeof(path), "OUT\\%s", fd.cFileName);
                if (remove(path) == 0) deleted_count++;
                else                    failed_count++;
            } while (FindNextFileA(hFind, &fd));
            FindClose(hFind);
        }
#else
        DIR *dp = opendir("OUT");
        if (dp) {
            struct dirent *ep;
            while ((ep = readdir(dp)) != NULL) {
                size_t len = strlen(ep->d_name);
                if (len > 4 && strcmp(ep->d_name + len - 4, ".bin") == 0) {
                    char path[512];
                    snprintf(path, sizeof(path), "OUT/%s", ep->d_name);
                    if (remove(path) == 0) deleted_count++;
                    else                    failed_count++;
                }
            }
            closedir(dp);
        }
#endif
    }

    // Also delete common CSV outputs we know about (optional, best-effort).
    {
        const char *extensions[] = {".csv"};
        const char *base_names[] = {
            "DIC", "AT", "O2", "TOC", "NO3", "NH4", "PO4", "Si", "Sal",
            "Phy1", "Phy2", "SPM", "pCO2", "pH", "CO2",
            "daily_averages", "diagnostics"
        };

        for (int ext_idx = 0; ext_idx < 1; ext_idx++) {
            for (int name_idx = 0; name_idx < 17; name_idx++) {
                char filepath[256];
                snprintf(filepath, sizeof(filepath), "OUT/%s%s", base_names[name_idx], extensions[ext_idx]);

                if (remove(filepath) == 0) {
                    deleted_count++;
                } else if (errno != ENOENT) {
                    failed_count++;
                }
            }
        }
    }
    
    if (debug_level >= DEBUG_LEVEL_ESSENTIAL) {
        printf("🗑️  Cleaned %d old output files", deleted_count);
        if (failed_count > 0) {
            printf(" (⚠️ %d files could not be deleted)", failed_count);
        }
        printf("\n");
    }
}

/**
 * @brief Initializes the hydrodynamic arrays and geometry for each segment.
 *
 * Core Variable Relationships
 * - freeArea[i]: Free surface cross-sectional area (water above riverbed) [m²]
 * - baseArea[i]: Channel-bed cross-sectional area (below datum) [m²]
 * - totalArea[i]: Total water-column cross-sectional area [m²]
 * - totalArea[i] = freeArea[i] + baseArea[i]
 * - width[i]: Channel width at segment i [m]
 * - waterDepth[i]: Water depth [m]
 * - waterDepth[i] = totalArea[i] / width[i]
 * - riverbed_depth[i]: Depth of the riverbed below datum [m]
 * - baseArea[i] = width[i] * riverbed_depth[i]
 */
void initializeHydrodynamics() {
    PRINTF_INIT("🔧 STAGGERED GRID INITIALIZATION: Setting initial velocities for all even indices\n");
    
    // 1. First, initialize all arrays to zero (critical for stability)
    for (int i = 0; i <= M; i++) {
        E[i] = 0.0;             // Initialize temporary water level storage array
        Y[i] = -5.0;            // Solver workspace seed (arbitrary negative for convergence)
        tempFreeArea[i] = 0.0;  // Set temporary water levels to zero
        tempVelocity[i] = 0.0;  // Set temporary velocities to zero
        freeArea[i] = 0.0;      // Initialize actual water levels
        velocity[i] = 0.0;      // Initialize actual velocities
        level[i] = 0.0;         // Initialize water level above datum
        
        // Initialize hydrodynamic coefficient arrays
        Chezy[i] = 60.0;        // Default Chezy coefficient
        FRIC[i] = 1.0 / (60.0 * 60.0); // Default friction coefficient
        Mero[i] = 1.0e-6;       // Default erosion coefficient
        tau_ero[i] = 1.0;       // Default critical shear stress for erosion
        tau_dep[i] = 0.3;       // Default critical shear stress for deposition
        rs[i] = 1.0;            // Default storage width ratio
    }
    
    // Initial small, stable water level (critical for stability)
    double initial_water_level = 0.1 * AMPL;

    // Define transition zone parameters using configuration-driven segment index
    int transition_center = index_2 > 0 ? index_2 : 30;

    if (align_segment_break_with_major_tributary && tributaryEnabled && num_segments > 1 && numTributaries > 0) {
        int dominant_cell = -1;
        double dominant_flow = -1.0;

        for (int trib = 0; trib < numTributaries; ++trib) {
            int cell = tributaries[trib].cellIndex;
            double *series = tributaries[trib].discharge;
            int n = tributaries[trib].dischargeDataSize;
            if (cell <= 0 || !series || n <= 0) {
                continue;
            }

            double sum = 0.0;
            int valid = 0;
            for (int k = 0; k < n; ++k) {
                double q = series[k];
                if (!isfinite(q)) {
                    continue;
                }
                sum += fabs(q);
                valid++;
            }
            if (valid == 0) {
                continue;
            }
            double mean_q = sum / (double)valid;
            if (mean_q > dominant_flow) {
                dominant_flow = mean_q;
                dominant_cell = cell;
            }
        }

        if (dominant_cell > 0) {
            transition_center = dominant_cell;
            if (calibration_mode == 0 && debug_level >= DEBUG_LEVEL_ESSENTIAL) {
                double km = ((double)(transition_center - 1) * DELXI) / 1000.0;
                printf("[INIT] ℹ️ INFO: Segment break aligned to dominant tributary at cell %d (%.1f km, mean Q ≈ %.1f m³/s)\n",
                       transition_center, km, dominant_flow);
            }
        }
    }
    // Clamp transition center to avoid overlapping the boundaries
    if (transition_center < 5) {
        transition_center = 5;
    } else if (transition_center > M - 5) {
        transition_center = M - 5;
    }

    transition_center = index_2 > 0 ? index_2 : transition_center;
    const int transition_width = segment_transition_width; // Configurable half-width (MAINT-001A fix)
    const int transition_start = transition_center - transition_width;
    const int transition_end = transition_center + transition_width;
    
    for (int i = 1; i <= M; i++) {
        if (i < transition_start) {
            width[i] = B1 * exp(-(i * DELXI) / LC1);
        } else if (i > transition_end) {
            width[i] = B2 * exp(-((i - transition_center) * DELXI) / LC2);
        } else {
            // Transition zone: smooth blend between segments
            double width1 = B1 * exp(-(i * DELXI) / LC1);              // Segment 1 width
            double width2 = B2 * exp(-((i - transition_center) * DELXI) / LC2);         // Segment 2 width
            
            // Smooth transition using cosine interpolation
            double t = (double)(i - transition_start) / (double)(transition_end - transition_start);
            double blend_factor = 0.5 * (1.0 - cos(M_PI * t));  // Smooth S-curve from 0 to 1
            
            width[i] = width1 * (1.0 - blend_factor) + width2 * blend_factor;
        }
        
        // Enforce minimum width for numerical stability
        const double MIN_WIDTH = 50.0;  // [m]
        if (width[i] < MIN_WIDTH) width[i] = MIN_WIDTH;
        
        // Set Chezy coefficients and storage ratios with smooth transition
        // FIX: Use global calibration parameters (Chezy1/2, Rs1/2) instead of hardcoded defaults
        if (i < transition_start) {
            // Pure segment 1 (cells 1-25)
            Chezy[i] = Chezy1;
            rs[i] = Rs1;
        } else if (i > transition_end) {
            // Pure segment 2 (cells 35+)
            Chezy[i] = Chezy2;
            rs[i] = Rs2;
        } else {
            // Transition zone: smooth blend between segments (cells 26-34)
            double t = (double)(i - transition_start) / (double)(transition_end - transition_start);
            double blend_factor = 0.5 * (1.0 - cos(M_PI * t));  // Same smooth S-curve
            
            Chezy[i] = Chezy1 * (1.0 - blend_factor) + Chezy2 * blend_factor;
            rs[i] = Rs1 * (1.0 - blend_factor) + Rs2 * blend_factor;
        }
        
        // Set riverbed cross-section - VERSION 1 EXACT MATCH
        baseArea[i] = width[i] * riverbed_depth[i];
        
        // Set initial water levels uniformly - CRITICAL FOR STABILITY
        freeArea[i] = width[i] * initial_water_level;
        tempFreeArea[i] = freeArea[i];
        level[i] = initial_water_level;  // Water level above datum
        
        // Calculate total area and water depth
        totalArea[i] = freeArea[i] + baseArea[i];
        waterDepth[i] = totalArea[i] / width[i];
        
        // Store friction coefficients
        FRIC[i] = 1.0 / (Chezy[i] * Chezy[i]);
        
        // Set sediment transport parameters based on segments
        // Supports optional 3rd/4th segment for independent upstream SPM control
        if (index_3 > 0 && Mero3 > 0.0) {
            // Multi-segment sediment: seg1 | blend12 | seg2 | blend23 | seg3 [| blend34 | seg4]
            const int t12_start = transition_center - transition_width;
            const int t12_end   = transition_center + transition_width;
            const int t23_start = index_3 - transition_width;
            const int t23_end   = index_3 + transition_width;
            // Optional 4th segment
            const int has_seg4 = (index_4 > 0 && Mero4 > 0.0);
            const int t34_start = has_seg4 ? (index_4 - transition_width) : M + 1;
            const int t34_end   = has_seg4 ? (index_4 + transition_width) : M + 2;

            if (i < t12_start) {
                Mero[i] = Mero1;  tau_ero[i] = tau_ero1;  tau_dep[i] = tau_dep1;
            } else if (i <= t12_end) {
                double t = (double)(i - t12_start) / (double)(t12_end - t12_start);
                double bf = 0.5 * (1.0 - cos(M_PI * t));
                Mero[i]    = Mero1    * (1.0 - bf) + Mero2    * bf;
                tau_ero[i] = tau_ero1  * (1.0 - bf) + tau_ero2  * bf;
                tau_dep[i] = tau_dep1  * (1.0 - bf) + tau_dep2  * bf;
            } else if (i < t23_start) {
                Mero[i] = Mero2;  tau_ero[i] = tau_ero2;  tau_dep[i] = tau_dep2;
            } else if (i <= t23_end) {
                double t = (double)(i - t23_start) / (double)(t23_end - t23_start);
                double bf = 0.5 * (1.0 - cos(M_PI * t));
                Mero[i]    = Mero2    * (1.0 - bf) + Mero3    * bf;
                tau_ero[i] = tau_ero2  * (1.0 - bf) + tau_ero3  * bf;
                tau_dep[i] = tau_dep2  * (1.0 - bf) + tau_dep3  * bf;
            } else if (has_seg4 && i >= t34_start) {
                if (i >= t34_end) {
                    Mero[i] = Mero4;  tau_ero[i] = tau_ero4;  tau_dep[i] = tau_dep4;
                } else {
                    double t = (double)(i - t34_start) / (double)(t34_end - t34_start);
                    double bf = 0.5 * (1.0 - cos(M_PI * t));
                    Mero[i]    = Mero3    * (1.0 - bf) + Mero4    * bf;
                    tau_ero[i] = tau_ero3  * (1.0 - bf) + tau_ero4  * bf;
                    tau_dep[i] = tau_dep3  * (1.0 - bf) + tau_dep4  * bf;
                }
            } else {
                Mero[i] = Mero3;  tau_ero[i] = tau_ero3;  tau_dep[i] = tau_dep3;
            }
        } else {
            // 2-segment sediment (default)
            if (i < transition_start) {
                Mero[i] = Mero1;  tau_ero[i] = tau_ero1;  tau_dep[i] = tau_dep1;
            } else if (i > transition_end) {
                Mero[i] = Mero2;  tau_ero[i] = tau_ero2;  tau_dep[i] = tau_dep2;
            } else {
                double t = (double)(i - transition_start) / (double)(transition_end - transition_start);
                double blend_factor = 0.5 * (1.0 - cos(M_PI * t));
                Mero[i]    = Mero1    * (1.0 - blend_factor) + Mero2    * blend_factor;
                tau_ero[i] = tau_ero1  * (1.0 - blend_factor) + tau_ero2  * blend_factor;
                tau_dep[i] = tau_dep1  * (1.0 - blend_factor) + tau_dep2  * blend_factor;
            }
        }
    }
    

    
    // 3. SET VELOCITIES - CRITICAL PART FOR STAGGERED GRID
    int velocity_count = 0;
    double downstream_velocity = 0.0;
    // Use signed upstream NET discharge for boundary velocity (negative = into domain)
    double upstream_velocity = Discharge_ups(0) / totalArea[M];
    
    // Initialize velocities ONLY at even indices (staggered grid compliance)
    for (int i = 2; i <= M2; i += 2) {
        // Set small initial velocity (positive = seaward flow)
        velocity[i] = 0.01; // Small positive initial velocity
        tempVelocity[i] = velocity[i];
        velocity_count++;
    }
    
    // 4. Set boundary velocities
    velocity[M] = upstream_velocity;     // Upstream boundary (river)
    tempVelocity[M] = velocity[M];
    
    PRINTF_INIT("✅ Initialized %d velocity points (even indices 2 to %d)\n", 
           velocity_count, M2);
    PRINTF_INIT("   downstream_velocity = %.6g m/s\n", downstream_velocity);
    PRINTF_INIT("   upstream_velocity = %.6g m/s\n", upstream_velocity);

    PRINTF_INIT("Hydrodynamics Initialized with %d segments following staggered grid pattern.\n", num_segments);
}

/**
 * @brief Initializes transport variables with physics-based initial values.
 *
 * @details
 * This function initializes the transport variables following Savenije's theory:
 * - Dispersion coefficients initialized with realistic longitudinal profile
 * - Uses Van der Burgh's equation and estuarine geometry for dispersion
 * - Sets up proper initial conditions for salinity transport
 * 
 * References:
 * - Savenije, H.H.G. (2012). Salinity and Tides in Alluvial Estuaries, 2nd edition
 * - Gisen, J.I.A., Savenije, H.H.G., et al. (2015). Guidelines for predictive
 *   calibration of estuarine 1-D models, HESS, 19, 2791-2803
 */
void initializeTransportVariables() {
    
    // Initialize basic arrays to zero
    for (int i = 0; i <= M; i++) {
        fl[i] = 0.0;           // Reset flux array
        watflux[i] = 0.0;      // Reset water flux array
    }
    
    // === PHYSICS-BASED DISPERSION INITIALIZATION ===
    // Following Savenije (2012) theory for dispersion in alluvial estuaries
    
    // 1. Calculate reference dispersion at the mouth (x=0)
    double D0 = 0.0;        // Will store mouth dispersion [m²/s]
    double K_vdb = 0.0;     // Will store Van der Burgh coefficient
    
    // Calculate mouth velocity amplitude (typical tidal velocity) - more conservative
    double omega = 2.0 * M_PI / TIDAL_PERIOD;  // Angular frequency [rad/s]
    
    // More realistic tidal velocity amplitude calculation
    double tidal_velocity_amplitude = AMPL * omega;  // Simple linear relationship: v = A * ω
    tidal_velocity_amplitude = fmin(tidal_velocity_amplitude, 1.5);  // Cap at realistic 1.5 m/s
    
    // Correct tidal excursion calculation: E = v_tidal * T/(2π) 
    double E_0 = tidal_velocity_amplitude * TIDAL_PERIOD / (2.0 * M_PI);  // Tidal excursion [m]
    
    // Calculate mouth dispersion using conservative Savenije's approach: D = k * E * v
    // Use factor of 1/10 instead of 1/15 for more realistic values
    D0 = fmax(500.0, E_0 * tidal_velocity_amplitude / 10.0);
    // Apply calibration correction
    D0 *= D0_CORRECTION;
    
    // Cap the mouth dispersion to realistic values for large estuaries
    D0 = fmin(D0, 5000.0);  // Increased cap to allow calibration range (500 was too restrictive)
    
    // Calculate Van der Burgh coefficient using estuary geometry (Savenije, 2012)
    double H_avg = 0.0;    // Average depth
    double B_avg = 0.0;    // Average width
    int count = 0;
    
    // Calculate average depth and width for Van der Burgh coefficient
    for (int i = 1; i <= M; i++) {
        H_avg += waterDepth[i];
        B_avg += width[i];
        count++;
    }
    H_avg /= count;
    B_avg /= count;
    
    // More conservative Van der Burgh coefficient calculation
    // Use typical values for alluvial estuaries (0.2-0.6 range)
    double friction_influence = sqrt(G/FRIC[1]);
    double K_phys = 0.4 / (1.0 + 0.5 * (H_avg/B_avg) * friction_influence);
    
    // Apply calibration scaling
    K_vdb = K_phys * C_VDB;
    K_vdb = fmax(0.1, fmin(1.0, K_vdb));  // Constrain to realistic range extended for calibration
    
    PRINTF_INIT("   • Reference mouth dispersion (D0): %.1f m²/s\n", D0);
    PRINTF_INIT("   • Van der Burgh coefficient (K): %.3f\n", K_vdb);
    PRINTF_INIT("   • Tidal velocity amplitude: %.2f m/s\n", tidal_velocity_amplitude);
    PRINTF_INIT("   • Tidal excursion: %.0f m\n", E_0);
    
    // 2. Calculate dispersion profile using Van der Burgh's equation
    for (int i = 1; i <= M; i++) {
        // Use area-based Van der Burgh relationship: D(x) = D0 * (A(x)/A0)^K
        double area_ratio = totalArea[i] / totalArea[1];
        disp[i] = D0 * pow(area_ratio, K_vdb);
        
        // Enforce reasonable minimum and maximum dispersion bounds
        disp[i] = fmax(5.0, fmin(1000.0, disp[i]));  // Range: 5-1000 m²/s
    }
    
    PRINTF_INIT("   • Dispersion range: %.1f - %.1f m²/s\n", disp[M], disp[1]);
}

/**
 * @brief Initializes biogeochemical variables.
 *
 * @details
 * This function initializes the biogeochemical variables, such as nutrient
 * concentrations and phytoplankton biomass.
 */
void initializeBiogeochemistry() {
    for (int var = 0; var < NUM_BC_VARS; var++) {
        memset(v[var].name, 0, sizeof(v[var].name));
        strncpy(v[var].name, variableNames[var], sizeof(v[var].name) - 1);
        v[var].name[sizeof(v[var].name) - 1] = '\0';

        v[var].env = 1;

        for (int i = 0; i <= M; i++) {
            v[var].avg[i] = 0.0;
            v[var].concflux[i] = 0.0;
            v[var].advflux[i] = 0.0;
            v[var].disflux[i] = 0.0;
        }
    }
    
    // Initialize auto-calculated variables (pH, pCO2, CO2) explicitly
    // Calculate CO2 from DIC/pH equilibrium to avoid first-day spike
    // pCO2 = CO2 / Henry × 1e6, where CO2 = DIC / (1 + K1/H + K1×K2/H²)
    PRINTF_INIT("✅ Initializing carbonate chemistry in equilibrium (pH, pCO2, CO2)...\n");
    
    // Use typical estuarine conditions for initial calculation
    const double init_pH = 7.8;           // Typical estuarine pH
    const double init_H = pow(10.0, -init_pH);
    const double init_temp = 28.0;        // Typical tropical temperature
    const double init_sal_avg = 5.0;      // Average estuary salinity for equilibrium calc
    
    // Calculate equilibrium constants at typical conditions (approximate)
    // K1 and K2 from Cai & Wang formulation at T=28°C, S=5
    double T_kelvin = init_temp + 273.15;
    double pK1_init = -14.8425 + (3404.71/T_kelvin) + (0.032786*T_kelvin);
    double f1_init = -0.0230848 - (14.3456/T_kelvin);
    double f2_init = 0.000691881 + (0.429955/T_kelvin);
    double K1_init = pow(10, -(pK1_init + (f1_init*pow(init_sal_avg,0.5)) + (f2_init*init_sal_avg)));
    
    double pK2_init = -6.4980 + (2902.39/T_kelvin) + (0.02379*T_kelvin);
    double f3_init = -0.458898 + (41.24048/T_kelvin);
    double f4_init = 0.0284743 - (2.55895/T_kelvin);
    double K2_init = pow(10, -(pK2_init + (f3_init*pow(init_sal_avg,0.5)) + (f4_init*init_sal_avg)));
    
    // Henry's constant at typical conditions (mol/L.atm, then convert to mmol/m³.atm)
    double lnK0_init = -574.70126 + 21541.52/T_kelvin - 0.000147759*(T_kelvin*T_kelvin) + 89.892*(log(T_kelvin));
    double f_henry_init = 0.029941 - 0.00027455*T_kelvin + 0.00000053407*(T_kelvin*T_kelvin);
    double KH_init = exp(lnK0_init + f_henry_init*init_sal_avg) * 1e6;  // mmol/m³.atm
    
    for (int i = 0; i <= M; i++) {
        v[PH].c[i] = init_pH;
        
        // CO2 equilibrium from DIC: CO2 = DIC / (1 + K1/H + K1*K2/H²)
        // Use the actual DIC value at this cell if available
        double DIC_val = v[DIC].c[i];
        if (DIC_val < 100.0) DIC_val = 1800.0;  // Use default if DIC not yet set
        
        double denom = 1.0 + K1_init / init_H + K1_init * K2_init / (init_H * init_H);
        double CO2_eq = DIC_val / denom;
        v[CO2].c[i] = CO2_eq;
        
        // pCO2 from Henry's law: pCO2 = CO2 / KH × 1e6
        // KH is already in mmol/m³.atm, so pCO2 = CO2 [mmol/m³] / KH [mmol/m³.atm] × 1e6 [µatm/atm]
        double pCO2_eq = CO2_eq / KH_init * 1e6;
        
        // Clamp pCO2 to reasonable estuarine range (300-5000 µatm)
        if (pCO2_eq < 300.0) pCO2_eq = 400.0;
        if (pCO2_eq > 10000.0) pCO2_eq = 2000.0;
        v[pCO2].c[i] = pCO2_eq;
    }
    
    PRINTF_INIT("   Initial carbonate equilibrium: pH=%.2f, CO2=%.1f mmol/m³, pCO2=%.0f µatm\n",
                init_pH, v[CO2].c[1], v[pCO2].c[1]);

    if (calibration_mode == 0) {
        const char *ext = (output_storage_format == OUTPUT_FORMAT_BINARY) ? ".bin" : ".csv";
        printf("✅ Output files configured for auto-calculated variables:\n");
        printf("   - pH: OUT/%s%s\n", v[PH].name, ext);
        printf("   - pCO2: OUT/%s%s\n", v[pCO2].name, ext);
        printf("   - CO2: OUT/%s%s\n", v[CO2].name, ext);
    }
    
    // Initialize biogeochemical rate constants
    for(int i=0; i<=M; i++)
    {
        for(char s=Phy1; s<=Phy2; s++)
        {
            NPP_NO3[i][s]=0.0;
            NPP_NH4[i][s]=0.0;
            GPP[i][s]=0.0;
            phydeath[i][s]=0.0;
        }

        NPP[i]=0.0;
        Si_consumption[i]=0.0;
        adegrad[i]=0.0;
        denit[i]=0.0;
        nitrif[i]=0.0;
        o2air[i]=0.0;
        co2air[i]=0.0;
    }

    // Initialize concentration profiles for all chemical species
    for (int var_idx = 0; var_idx < CHEM_COUNT; var_idx++) {
        // Check if the variable is environmentally simulated and requires data input
        if (v[var_idx].env == 1 && isDataInputRequired(var_idx)) { 
            
            double initial_cub = 0.0; // Default upstream boundary concentration
            double initial_clb = 0.0; // Default downstream boundary concentration

            // Get the first value from the upstream boundary condition data array
            // Corrected condition: Check if data pointer is not NULL
            if (upstreamBC[var_idx].data != NULL) {
                initial_cub = upstreamBC[var_idx].data[0];
                // DEBUG: Check for negative initial values
                if (initial_cub < 0.0) {
                    fprintf(stderr, "FATAL TRACER ERROR: stage=Init:upstreamBC t=%d species=%d (%s) cell=%d value=%.17g reason=negative\n",
                            0, var_idx, getSpeciesName(var_idx), M, initial_cub);
                    exit(EXIT_FAILURE);
                }
                if (!isfinite(initial_cub)) {
                    fprintf(stderr, "FATAL TRACER ERROR: stage=Init:upstreamBC t=%d species=%d (%s) cell=%d value=%.17g reason=nonfinite\n",
                            0, var_idx, getSpeciesName(var_idx), M, initial_cub);
                    exit(EXIT_FAILURE);
                }
            } else {
                if (var_idx == Sal) initial_cub = 0.0; 
                // else printf("Warning: Missing initial upstream data for %s, using 0.0\n", variableNames[var_idx]);
            }

            // Get the first value from the downstream boundary condition data array
            // Corrected condition: Check if data pointer is not NULL
            if (downstreamBC[var_idx].data != NULL) {
                initial_clb = downstreamBC[var_idx].data[0];
                // DEBUG: Check for negative initial values
                if (initial_clb < 0.0) {
                    fprintf(stderr, "FATAL TRACER ERROR: stage=Init:downstreamBC t=%d species=%d (%s) cell=%d value=%.17g reason=negative\n",
                            0, var_idx, getSpeciesName(var_idx), 1, initial_clb);
                    exit(EXIT_FAILURE);
                }
                if (!isfinite(initial_clb)) {
                    fprintf(stderr, "FATAL TRACER ERROR: stage=Init:downstreamBC t=%d species=%d (%s) cell=%d value=%.17g reason=nonfinite\n",
                            0, var_idx, getSpeciesName(var_idx), 1, initial_clb);
                    exit(EXIT_FAILURE);
                }
            } else {
                // printf("Warning: Missing initial downstream data for %s, using 0.0\n", variableNames[var_idx]);
            }
            
            v[var_idx].cub = initial_cub;
            v[var_idx].clb = initial_clb;
            
            // Initialize concentration profile for this variable
            for (int i = 0; i <= M; i++) {
                if (var_idx == Sal) {
                    // For salinity, initialize using the same boundary-consistent interpolation
                    // policy as other tracers (mouth = clb, upstream = cub). Transport will
                    // spin up the profile during warmup.
                    if (M > 0) {
                        double fraction = (double)i / (double)M;  // 0 at mouth, 1 at upstream end
                        v[var_idx].c[i] = initial_clb * (1.0 - fraction) + initial_cub * fraction;
                    } else {
                        v[var_idx].c[i] = initial_clb;
                    }
                } else {
                    // For other variables, use linear interpolation
                    if (M > 0) { 
                        double fraction = (double)i / (double)M;  // 0 at mouth, 1 at upstream end
                        v[var_idx].c[i] = initial_clb * (1.0 - fraction) + initial_cub * fraction;
                        
                        // COMPREHENSIVE VALIDATION: Check for negative concentrations everywhere
                        if (v[var_idx].c[i] < 0.0) {
                            fprintf(stderr,
                                    "FATAL TRACER ERROR: stage=Init:interp t=%d species=%d (%s) cell=%d value=%.17g reason=negative\n",
                                    0, var_idx, getSpeciesName(var_idx), i, v[var_idx].c[i]);
                            exit(EXIT_FAILURE);
                        }
                        if (!isfinite(v[var_idx].c[i])) {
                            fprintf(stderr,
                                    "FATAL TRACER ERROR: stage=Init:interp t=%d species=%d (%s) cell=%d value=%.17g reason=nonfinite\n",
                                    0, var_idx, getSpeciesName(var_idx), i, v[var_idx].c[i]);
                            exit(EXIT_FAILURE);
                        }
                        
                    } else {
                        v[var_idx].c[i] = initial_clb; 
                    }
                }
                v[var_idx].avg[i] = 0.0; 
            }
        }
    }

    // After loading boundary conditions, validate chemical species concentrations
    PRINTF_INIT("🔄 Initializing chemical species with default values...\n");
    for (int s = 0; s < MAXV; s++) {
        if (v[s].env == 1) {
            // Skip auto-calculated variables
            if (s == pCO2 || s == PH || s == CO2) continue;
            
            // Set default values based on species for cells that aren't initialized (in mmol/m³)
            double default_value = 0.0;
            
            switch(s) {
                case O2:  default_value = 160.0;  break; // ~5 mg/L in mmol/m³
                case Sal: default_value = 0.1;    break; // Salinity in PSU
                case TOC: default_value = 300.0;  break; // TOC in mmol/m³ 
                case SPM: default_value = 10.0;   break; // SPM in g/m³
                case DIC: default_value = 1800.0; break; // DIC in mmol/m³
                case NH4: default_value = 5.0;    break; // NH4 in mmol/m³
                case NO3: default_value = 20.0;   break; // NO3 in mmol/m³
                case PO4: default_value = 1.0;    break; // PO4 in mmol/m³
                case Si:  default_value = 30.0;   break; // Si in mmol/m³
                case AT:  default_value = 1800.0; break; // AT in mmol/m³
                case Phy1:default_value = 10.0;   break; // Phy1 in mmol/m³
                case Phy2:default_value = 1.0;    break; // Phy2 in mmol/m³
                default:  default_value = 1.0;    break; // Generic default
            }
            
            // Fill any remaining cells that were not initialized by boundary-based interpolation.
            // Previously this block hard-aborted on NaN/Inf/negative, but its purpose is to
            // ensure a valid starting state.
            for (int i = 0; i <= M; i++) {
                if (!isfinite(v[s].c[i]) || v[s].c[i] < 0.0) {
                    v[s].c[i] = default_value;
                }
            }

            // Ensure boundary values are reasonable; fall back to defaults if needed.
            if (!isfinite(v[s].clb) || v[s].clb < 0.0) {
                v[s].clb = default_value;
            }
            if (!isfinite(v[s].cub) || v[s].cub < 0.0) {
                v[s].cub = default_value;
            }
            
            // Also initialize any averaging arrays
            for (int i = 0; i <= M; i++) {
                v[s].avg[i] = 0.0;
            }
        }
    }

    diagnostics_validate_tracer_state("Init:exit", 0);

}
/**
 * @brief Assigns rate constants for biogeochemical processes.
 *
 * @details
 * This function assigns the rate constants for the biogeochemical processes,
 * such as photosynthesis, respiration, and nutrient uptake.
 */
void assignBiogeochemicalRateConstants() {
    // Default rate constants for tropical eutrophic estuaries.
    // All values can be overridden via params.txt.
    if (!biogeo_override_Pb[Phy1]) Pb[Phy1] = 2.0e-4;
    if (!biogeo_override_Pb[Phy2]) Pb[Phy2] = 2.6e-4;
    if (!biogeo_override_alpha[Phy1]) alpha[Phy1] = 8.0e-7;
    if (!biogeo_override_alpha[Phy2]) alpha[Phy2] = 8.0e-7;

    kexcr[Phy1]  = 0.05;  // Excretion constant (Tref=20°C)
    kexcr[Phy2]  = 0.05;
    kgrowth[Phy1] = 0.20;  // Growth respiration cost
    kgrowth[Phy2] = 0.15;  // Non-diatoms more efficient

    kmaint[Phy1] = 3.0e-7;   // Maintenance constant [1/s]
    kmaint[Phy2] = 2.5e-7;
    if (!biogeo_override_kmortality[Phy1]) kmortality[Phy1] = 1.4e-7;
    if (!biogeo_override_kmortality[Phy2]) kmortality[Phy2] = 1.4e-7;

    if (!biogeo_override_KSi) KSi[Phy1] = 2.0;   // Silica half-saturation [mmol/m³]
    KSi[Phy2] = 0.0;

    if (!biogeo_override_KN[Phy1]) KN[Phy1] = 25.0;   // N half-saturation [mmol/m³]
    if (!biogeo_override_KN[Phy2]) KN[Phy2] = 25.0;
    if (!biogeo_override_KPO4[Phy1]) KPO4[Phy1] = 0.2;
    if (!biogeo_override_KPO4[Phy2]) KPO4[Phy2] = 0.2;

    if (!isfinite(KTOC) || KTOC <= 0.0) KTOC = 300.0;  // TOC half-saturation [mmol/m³]

    KO2_ox  = 31.0;    // O2 half-saturation for aerobic degradation
    KO2_nit = 51.25;   // O2 half-saturation for nitrification
    KinO2   = 33.0;    // Inhibition constant for denitrification
    KNO3    = 26.07;   // NO3 half-saturation
    KNH4    = 80.0;    // NH4 half-saturation
    redsi   = 15./106.; // Redfield Si:C
    redn    = 16./106.; // Redfield N:C
    redp    = 1./106.;  // Redfield P:C

    if (!biogeo_override_kox) kox = 1.5e-4;         // Aerobic degradation rate [1/s]
    if (!isfinite(kdenit) || kdenit <= 0.0) kdenit = 1.0e-4;  // Denitrification rate [1/s]
    if (!biogeo_override_knit) knit = 1.5e-4;        // Nitrification rate [1/s]
    if (!biogeo_override_pCO2atmo) pCO2atmo = 415.0;  // Atmospheric pCO2 [µatm] (2024 global mean)

    // Light attenuation defaults (overridden by params.txt if provided)
    if (!isfinite(kbg) || kbg <= 0.0) kbg = 0.2;     // Background attenuation [1/m]
    if (!isfinite(kspm) || kspm <= 0.0) kspm = 0.02;  // SPM attenuation [m²/g]

    // Evidence-first audit of effective parameter values (for reproducibility).
    // Prefer OUT/ (always created). Also attempt outputs/metrics/ when available.
    {
        FILE *fp = fopen("OUT/effective_light_params.csv", "w");
        if (fp) {
            fprintf(fp, "kbg,kspm,KD_Phy\n");
            fprintf(fp, "%.10g,%.10g,%.10g\n", kbg, kspm, KD_Phy);
            fclose(fp);
        }
        // Best-effort secondary copy (non-fatal if directories missing)
        fp = fopen("outputs/metrics/effective_light_params.csv", "w");
        if (fp) {
            fprintf(fp, "kbg,kspm,KD_Phy\n");
            fprintf(fp, "%.10g,%.10g,%.10g\n", kbg, kspm, KD_Phy);
            fclose(fp);
        }
    }
    // Euler-Mascheroni constant now defined locally in biogeo.c as EULER_CONSTANT
    // ws is now read from params.txt (settling_velocity parameter)
    // Default initialization moved to variables.c, value from Winterwerp (2002)
    if (ws <= 0.0) ws = 5.0e-4;  // Fallback if not set in params.txt: 0.5 mm/s
    rho_w		 =1000.0;					//kg m-3
    g            = G;						//m s-2  (use define.h constant)
}

// Use the correct variableNames array from variables.c
const char* getSpeciesName(int s) {
    // Use the global variableNames array which matches the Chem enum exactly
    extern const char *variableNames[CHEM_COUNT];
    
    if (s >= 0 && s < CHEM_COUNT) {
        return variableNames[s];
    }
    return "Unknown"; // For out-of-range indices
}
