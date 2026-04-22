#include "file.h"
#include "define.h"
#include "variables.h"
#include "diagnostics.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <errno.h>

/**
 * @brief Lightweight parser to retrieve calibration flags before full initialization.
 */
int read_calibration_flags(const char *filename,
                           int *calibration_mode_out,
                           int *debug_level_out) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        return -1;
    }

    int mode_local = calibration_mode_out ? *calibration_mode_out : 0;
    int debug_local = debug_level_out ? *debug_level_out : debug_level;
    bool mode_found = false;
    bool debug_found = false;

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '#' || line[0] == '\n') {
            continue;
        }

        if (!mode_found) {
            int parsed_mode = 0;
            if (sscanf(line, "calibration_mode = %d", &parsed_mode) == 1) {
                if (parsed_mode < 0) parsed_mode = 0;
                if (parsed_mode > 7) parsed_mode = 7;
                mode_local = parsed_mode;
                mode_found = true;
                continue;
            }
        }

        if (!debug_found) {
            int parsed_debug = 0;
            if (sscanf(line, "debug_level = %d", &parsed_debug) == 1) {
                debug_local = parsed_debug;
                debug_found = true;
                continue;
            }
        }

        if (mode_found && debug_found) {
            break;
        }
    }

    fclose(file);

    // Strict scientific configuration: no environment-variable overrides for runtime clocks.

    if (calibration_mode_out) {
        *calibration_mode_out = mode_local;
    }
    if (debug_level_out) {
        *debug_level_out = debug_local;
    }
    return 0;
}

// -----------------------------------------------------------------------------
// Chemical boundary mapping helpers
// -----------------------------------------------------------------------------

typedef struct {
    double *upstreamTime;
    double *upstreamData;
    double *downstreamTime;
    double *downstreamData;
    int required;
} ChemicalBoundaryData;

static ChemicalBoundaryData chemicalBoundaryMap[NUM_BC_VARS];

static void setupChemicalBoundaryDataMap(void) {
    for (int i = 0; i < CHEM_COUNT; i++) {
        chemicalBoundaryMap[i].upstreamTime = upstreamBC[i].time;
        chemicalBoundaryMap[i].upstreamData = upstreamBC[i].data;
        chemicalBoundaryMap[i].downstreamTime = downstreamBC[i].time;
        chemicalBoundaryMap[i].downstreamData = downstreamBC[i].data;
        chemicalBoundaryMap[i].required = 1;
    }
}

// -----------------------------------------------------------------------------
// Runtime overrides (configuration-driven, no INPUT file edits required)
// -----------------------------------------------------------------------------

/**
 * @brief Sanitize a tributary name to an env-var-friendly identifier.
 *
 * Converts to uppercase and replaces non-alphanumeric characters with '_'.
 * Example: "MyTrib" -> "MYTRIB".
 */
static void sanitize_env_ident(const char *name, char *out, size_t out_cap)
{
    if (!out || out_cap == 0) {
        return;
    }
    out[0] = '\0';
    if (!name) {
        return;
    }

    size_t k = 0;
    for (const unsigned char *p = (const unsigned char *)name; *p && k + 1 < out_cap; ++p) {
        unsigned char c = *p;
        if (isalnum(c)) {
            out[k++] = (char)toupper(c);
        } else {
            // Collapse repeated underscores
            if (k == 0 || out[k - 1] == '_') {
                continue;
            }
            out[k++] = '_';
        }
    }
    // Trim trailing underscore
    if (k > 0 && out[k - 1] == '_') {
        k--;
    }
    out[k] = '\0';
}

/**
 * @brief Apply optional env-var overrides for tributary locations.
 *
 * Supported:
 *   - CGEM_TRIB_<NAME>_CELLINDEX (1..M)
 * Example:
 *   CGEM_TRIB_MYTRIB_CELLINDEX=66
 */
static void apply_tributary_location_overrides(void)
{
    if (!tributaries || numTributaries <= 0 || M <= 0) {
        return;
    }

    for (int j = 0; j < numTributaries; ++j) {
        char ident[64];
        sanitize_env_ident(tributaries[j].name, ident, sizeof(ident));
        if (ident[0] == '\0') {
            continue;
        }

        char env_key[128];
        snprintf(env_key, sizeof(env_key), "CGEM_TRIB_%s_CELLINDEX", ident);
        const char *raw = getenv(env_key);
        if (!raw || !*raw) {
            continue;
        }

        char *endp = NULL;
        errno = 0;
        long v = strtol(raw, &endp, 10);
        if (errno != 0 || endp == raw || (*endp && !isspace((unsigned char)*endp))) {
            printf("⚠️ Ignoring %s='%s' (expected integer cell index)\n", env_key, raw);
            continue;
        }

        if (v < 1 || v > (long)M) {
            printf("⚠️ Ignoring %s='%s' (valid range: 1..%d)\n", env_key, raw, M);
            continue;
        }

        int old = tributaries[j].cellIndex;
        tributaries[j].cellIndex = (int)v;
        if (debug_level >= DEBUG_LEVEL_ESSENTIAL) {
            printf("[INIT] ℹ️ INFO: %s location override: cellIndex %d -> %d\n",
                   tributaries[j].name, old, tributaries[j].cellIndex);
        }
    }
}

// -----------------------------------------------------------------------------
// Generic two-column reader and path validation helpers
// -----------------------------------------------------------------------------

static void trim_line_inplace(char *line) {
    if (!line) {
        return;
    }

    // Trim leading spaces/tabs
    char *start = line;
    while (*start == ' ' || *start == '\t') {
        start++;
    }
    if (start != line) {
        memmove(line, start, strlen(start) + 1);
    }

    // Trim trailing whitespace + CR/LF
    size_t n = strlen(line);
    while (n > 0) {
        char c = line[n - 1];
        if (c == '\n' || c == '\r' || c == ' ' || c == '\t') {
            line[n - 1] = '\0';
            n--;
            continue;
        }
        break;
    }
}

int readFile(int datamax, double *gg, double *ff, const char *s) {
    if (s == NULL || !*s) {
        printf("❌ Error: File path is empty.\n");
        return 0;
    }

    FILE *fptr1 = fopen(s, "r");
    if (!fptr1) {
        int err = errno;
        fprintf(stderr, "❌ Error opening file '%s': %s\n", s, strerror(err));
        return 0;
    }

    double data1out, data2out;
    char buffer[256];
    int i = 0, validDataPoints = 0;
    int nearZeroCount = 0, negativeCount = 0, jumpCount = 0;
    double totalValue = 0.0, previousValue = 0.0;
    const double jumpThreshold = 10.0;

    int format_type = 0; // 0 = unknown, 1 = comma-separated, 2 = space/tab-separated
    bool has_header = false;

    const int isLightFile = (forcingData[FORCING_LIGHT].filePath &&
                             strcmp(s, forcingData[FORCING_LIGHT].filePath) == 0);

    if (fgets(buffer, sizeof(buffer), fptr1)) {
        if (buffer[0] == '#' ||
            buffer[0] == '%' ||
            strstr(buffer, "Time") ||
            strstr(buffer, "time") ||
            strstr(buffer, "DATE")) {
            has_header = true;
        } else {
            if (strchr(buffer, ',')) {
                format_type = 1;
                if (sscanf(buffer, "%lf,%lf", &data1out, &data2out) == 2) {
                    gg[i] = data1out;
                    ff[i] = data2out;
                    totalValue += data2out;
                    validDataPoints++;
                    previousValue = data2out;
                    i++;
                } else {
                    printf("⚠️ Warning: Could not parse first line as CSV in %s: %s", s, buffer);
                }
            } else {
                format_type = 2;
                if (sscanf(buffer, "%lf %lf", &data1out, &data2out) == 2) {
                    gg[i] = data1out;
                    ff[i] = data2out;
                    totalValue += data2out;
                    validDataPoints++;
                    previousValue = data2out;
                    i++;
                } else {
                    printf("⚠️ Warning: Could not parse first line as space-separated in %s: %s", s, buffer);
                }
            }
        }
    }

    if (has_header && format_type == 0) {
        if (fgets(buffer, sizeof(buffer), fptr1)) {
            format_type = strchr(buffer, ',') ? 1 : 2;
        }
    }

    if (format_type == 0) {
        format_type = 1;
    }

    while (i < datamax && fgets(buffer, sizeof(buffer), fptr1)) {
        if (buffer[0] == '#' || buffer[0] == '%' || buffer[0] == '\n') {
            continue;
        }

        bool parsed = false;
        if (format_type == 1) {
            if (sscanf(buffer, "%lf,%lf", &data1out, &data2out) == 2) {
                parsed = true;
            }
        } else {
            if (sscanf(buffer, "%lf %lf", &data1out, &data2out) == 2) {
                parsed = true;
            }
        }

        if (!parsed) {
            if (format_type == 1) {
                if (sscanf(buffer, "%lf %lf", &data1out, &data2out) == 2) {
                    parsed = true;
                    format_type = 2;
                    printf("⚠️ Switched to space-separated format mid-file in %s\n", s);
                }
            } else {
                if (sscanf(buffer, "%lf,%lf", &data1out, &data2out) == 2) {
                    parsed = true;
                    format_type = 1;
                    printf("⚠️ Switched to comma-separated format mid-file in %s\n", s);
                }
            }

            if (!parsed) {
                printf("⚠️ Skipping invalid line %d in %s: %s", i + 1, s, buffer);
                continue;
            }
        }

        totalValue += data2out;
        validDataPoints++;

        gg[i] = data1out;
        ff[i] = data2out;

        if (i > 0) {
            if (!isLightFile && data2out != 0 && fabs(previousValue) > 0.0) {
                double changeRatio = fabs(data2out / previousValue);
                if (changeRatio > jumpThreshold) {
                    jumpCount++;
                }
            }
        }

        if (data2out < 0) {
            negativeCount++;
        }

        previousValue = data2out;
        i++;
    }

    double averageValue = (validDataPoints > 0) ? (totalValue / validDataPoints) : 0.0;
    double nearZeroThreshold = 0.05 * fabs(averageValue);

    for (int j = 0; j < i; j++) {
        if (fabs(ff[j]) < nearZeroThreshold) {
            nearZeroCount++;
        }
    }

    fclose(fptr1);

    if (nearZeroCount > 0 || negativeCount > 0 || jumpCount > 0) {
        PRINTF_DATA("📊 Attention %s: %d zero (threshold: %.2f), %d negative, %d large jumps\n",
                    s, nearZeroCount, nearZeroThreshold, negativeCount, jumpCount);
    }

    return i;
}

int fileExists(const char *path) {
    FILE *file = fopen(path, "r");
    if (file) {
        fclose(file);
        return 1;
    }
    return 0;
}

int isValidPath(const char *path) {
    return (path && strlen(path) > 0);
}

void checkFilePath(const char *filePath, const char *varName) {
    if (!isValidPath(filePath)) {
        printf("❌ Warning: %s path is empty!\n", varName);
    } else if (!fileExists(filePath)) {
        printf("❌ Warning: %s does not exist: %s\n", varName, filePath);
    }
}

// -----------------------------------------------------------------------------
// Parameter and configuration parsing
// -----------------------------------------------------------------------------

void read_parameters(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("❌ Error: Cannot open parameter file");
        exit(EXIT_FAILURE);
    }

    // Optional crash diagnostics: stderr markers so we can localize faults
    // even when stdout is suppressed (e.g. during calibration).
    const char *trace_env = getenv("CGEM_TRACE_READ_PARAMS");
    const int trace_read_params = (trace_env && trace_env[0] != '\0' && atoi(trace_env) != 0);
    if (trace_read_params) {
        fprintf(stderr, "[TRACE] read_parameters: opened %s\n", filename ? filename : "(null)");
        fflush(stderr);
    }

    char line[256];
    num_segments = 0;

    // Phase 4/5: reset override flags and policy so repeated read_parameters() calls are deterministic.
    for (int s = 0; s < 2; ++s) {
        biogeo_override_Pb[s] = 0;
        biogeo_override_alpha[s] = 0;
        biogeo_override_kmortality[s] = 0;
        biogeo_override_KN[s] = 0;
        biogeo_override_KPO4[s] = 0;
    }
    biogeo_override_kox = 0;
    biogeo_override_knit = 0;
    biogeo_override_pCO2atmo = 0;
    biogeo_override_KSi = 0;
    carbonate_solver_failure_policy = CARBONATE_FAIL_HOLD_LAST;

    PRINTF_DATA("🔍 Reading parameters from %s...\n", filename);

    long dbg_line_no = 0;
    while (fgets(line, sizeof(line), file)) {
        dbg_line_no++;
        if (line[0] == '#' || line[0] == '\n') {
            continue;
        }

        if (trace_read_params && dbg_line_no <= 8) {
            // Show a few first non-empty lines to confirm encoding/format.
            fprintf(stderr, "[TRACE] read_parameters: line%ld='%.*s'\n", dbg_line_no, 160, line);
            fflush(stderr);
        }

        int parsed_mode = 0;
        if (sscanf(line, "calibration_mode = %d", &parsed_mode) == 1) {
            if (parsed_mode < 0) {
                parsed_mode = 0;
            } else if (parsed_mode > 7) {
                parsed_mode = 7;
            }
            calibration_mode = parsed_mode;
        }

        int parsed_stage = 0;
        if (sscanf(line, "calibration_stage = %d", &parsed_stage) == 1) {
            if (parsed_stage < 1) parsed_stage = 1;
            if (parsed_stage > 4) parsed_stage = 4;
            calibration_stage = parsed_stage;
        }

        char output_format_buf[32];
        if (sscanf(line, "output_storage_format = %31s", output_format_buf) == 1) {
            for (char *c = output_format_buf; *c; ++c) {
                *c = (char)toupper((unsigned char)*c);
            }
            if (strcmp(output_format_buf, "BINARY") == 0) {
                output_storage_format = OUTPUT_FORMAT_BINARY;
            } else if (strcmp(output_format_buf, "CSV") == 0) {
                output_storage_format = OUTPUT_FORMAT_CSV;
            } else {
                fprintf(stderr, "Warning: Unknown output_storage_format '%s'. Defaulting to CSV.\n", output_format_buf);
                output_storage_format = OUTPUT_FORMAT_CSV;
            }
        }

        int tmp_output_flag = 0;
        if (sscanf(line, "enable_flux_output = %d", &tmp_output_flag) == 1) {
            enable_flux_output = (tmp_output_flag != 0);
        }
        if (sscanf(line, "enable_reaction_output = %d", &tmp_output_flag) == 1) {
            enable_reaction_output = (tmp_output_flag != 0);
        }

        double calibration_tmp_days = 0.0;
        if (sscanf(line, "calibration_warmup_days = %lf", &calibration_tmp_days) == 1) {
            if (!isfinite(calibration_tmp_days) || calibration_tmp_days < 1.0) {
                printf("Warning: calibration_warmup_days must be >= 1 day; retaining %.1f d\n", calibration_warmup_days);
            } else {
                calibration_warmup_days = calibration_tmp_days;
            }
        }
        if (sscanf(line, "calibration_diagnostic_days = %lf", &calibration_tmp_days) == 1) {
            if (!isfinite(calibration_tmp_days) || calibration_tmp_days < 1.0) {
                printf("Warning: calibration_diagnostic_days must be >= 1 day; retaining %.1f d\n", calibration_diagnostic_days);
            } else {
                calibration_diagnostic_days = calibration_tmp_days;
            }
        }

        int calibration_tmp_iter = 0;
        if (sscanf(line, "calibration_max_iterations = %d", &calibration_tmp_iter) == 1) {
            if (calibration_tmp_iter < 1) {
                printf("Warning: calibration_max_iterations must be >= 1; retaining %d\n", calibration_max_iterations);
            } else {
                calibration_max_iterations = calibration_tmp_iter;
            }
        }


        // Warmup checkpointing settings removed (model always recomputes warmup).

        int temp_debug_level;
        if (sscanf(line, "debug_level = %d", &temp_debug_level) == 1) {
            if (temp_debug_level >= DEBUG_LEVEL_OFF && temp_debug_level <= DEBUG_LEVEL_DETAIL) {
                debug_level = temp_debug_level;
            } else {
                fprintf(stderr, "Warning: Invalid debug_level %d in params.txt, using default level 1 (ESSENTIAL)\n", temp_debug_level);
                debug_level = DEBUG_LEVEL_ESSENTIAL;
            }
        }

        if (sscanf(line, "MAXT = %ld", &MAXT) == 1) {}
        if (sscanf(line, "WARMUP = %ld", &WARMUP) == 1) {}
        if (sscanf(line, "DELTI = %d", &DELTI) == 1) {}
        if (sscanf(line, "TS = %d", &TS) == 1) {}
        if (sscanf(line, "DELXI = %d", &DELXI) == 1) {}
        if (sscanf(line, "EL = %d", &EL) == 1) {}
        if (sscanf(line, "AMPL = %lf", &AMPL) == 1) {}
        if (sscanf(line, "AMPL_CORRECTION = %lf", &AMPL_CORRECTION) == 1) {}

        int tmp_flag = 0;
        if (sscanf(line, "enable_upstream_discharge_smoothing = %d", &tmp_flag) == 1) {
            enable_upstream_discharge_smoothing = (tmp_flag != 0);
        }
        double tmp_tau_days = 0.0;
        if (sscanf(line, "upstream_discharge_smoothing_timescale_days = %lf", &tmp_tau_days) == 1) {
            if (!isfinite(tmp_tau_days) || tmp_tau_days < 0.0) {
                printf("Warning: upstream_discharge_smoothing_timescale_days must be >= 0; retaining %.2f d\n",
                       upstream_discharge_smoothing_timescale_days);
            } else {
                upstream_discharge_smoothing_timescale_days = tmp_tau_days;
            }
        }

        double tmp_discharge_scale = 0.0;
        if (sscanf(line, "upstream_discharge_scale = %lf", &tmp_discharge_scale) == 1) {
            if (!isfinite(tmp_discharge_scale) || tmp_discharge_scale <= 0.0) {
                printf("Warning: upstream_discharge_scale must be > 0; retaining %.3g\n", upstream_discharge_scale);
            } else {
                upstream_discharge_scale = tmp_discharge_scale;
            }
        }

        if (sscanf(line, "enable_tributary_subdaily_pulses = %d", &tmp_flag) == 1) {
            enable_tributary_subdaily_pulses = (tmp_flag != 0);
        }
        if (sscanf(line, "tributary_pulse_shape = %lf", &tmp_tau_days) == 1) {
            if (isfinite(tmp_tau_days) && tmp_tau_days >= 0.0) {
                tributary_pulse_shape = tmp_tau_days;
            }
        }
        if (sscanf(line, "urban_tributary_pulse_boost = %lf", &tmp_tau_days) == 1) {
            if (isfinite(tmp_tau_days) && tmp_tau_days > 0.0) {
                urban_tributary_pulse_boost = tmp_tau_days;
            }
        }

        char tmp_injection_mode[64];
        if (sscanf(line, "tributary_mass_injection_mode = %63s", tmp_injection_mode) == 1) {
            for (char *c = tmp_injection_mode; *c; ++c) {
                *c = (char)toupper((unsigned char)*c);
            }
            if (strcmp(tmp_injection_mode, "FISCHER_MASS_SOURCE") == 0) {
                tributary_mass_injection_mode = TRIB_INJECTION_FISCHER;
            } else if (strcmp(tmp_injection_mode, "RUTHERFORD_SPREAD") == 0) {
                tributary_mass_injection_mode = TRIB_INJECTION_RUTHERFORD;
            } else {
                fprintf(stderr,
                        "Warning: Unknown tributary_mass_injection_mode '%s'; using FISCHER_MASS_SOURCE.\n",
                        tmp_injection_mode);
                tributary_mass_injection_mode = TRIB_INJECTION_FISCHER;
            }
        }

        if (sscanf(line, "enable_biogeochemical_reactions = %d", &enable_biogeochemical_reactions) == 1) {}
        if (sscanf(line, "enable_carbonate_diagnostics = %d", &enable_carbonate_diagnostics) == 1) {}
        if (sscanf(line, "enable_suspended_sediment_dynamics = %d", &enable_suspended_sediment_dynamics) == 1) {}

        // Quadratic mortality closure (implicit zooplankton grazing)
        // Fennel et al. (2006, JGR-Oceans), Edwards & Yool (2000)
        double tmp_val_mort = 0.0;
        if (sscanf(line, "kmort2_Phy1 = %lf", &tmp_val_mort) == 1) {
            if (isfinite(tmp_val_mort) && tmp_val_mort >= 0.0) kmort2_Phy1 = tmp_val_mort;
        }
        // Benthic grazing clearance rate (Cloern 1982; Alpine & Cloern 1992)
        if (sscanf(line, "k_benth_graze = %lf", &tmp_val_mort) == 1) {
            if (isfinite(tmp_val_mort) && tmp_val_mort >= 0.0) k_benth_graze = tmp_val_mort;
        }
        if (sscanf(line, "kmort2_Phy2 = %lf", &tmp_val_mort) == 1) {
            if (isfinite(tmp_val_mort) && tmp_val_mort >= 0.0) kmort2_Phy2 = tmp_val_mort;
        }

        // Salinity stress mortality (Garnier et al., RIVE; Billen et al., 2001)
        if (sscanf(line, "kmort_sal = %lf", &tmp_val_mort) == 1) {
            if (isfinite(tmp_val_mort) && tmp_val_mort >= 0.0) kmort_sal = tmp_val_mort;
        }
        if (sscanf(line, "sal_tol_Phy1 = %lf", &tmp_val_mort) == 1) {
            if (isfinite(tmp_val_mort) && tmp_val_mort > 0.0) sal_tol_Phy1 = tmp_val_mort;
        }
        if (sscanf(line, "sal_tol_Phy2 = %lf", &tmp_val_mort) == 1) {
            if (isfinite(tmp_val_mort) && tmp_val_mort > 0.0) sal_tol_Phy2 = tmp_val_mort;
        }

        // ---------------------------------------------------------------------
        // Phase 4: Biogeochemical / carbonate parameter overrides
        // ---------------------------------------------------------------------
        double tmp_val = 0.0;
        if (sscanf(line, "Pb_Phy1 = %lf", &tmp_val) == 1) {
            Pb[Phy1] = tmp_val;
            biogeo_override_Pb[Phy1] = 1;
        }
        if (sscanf(line, "Pb_Phy2 = %lf", &tmp_val) == 1) {
            Pb[Phy2] = tmp_val;
            biogeo_override_Pb[Phy2] = 1;
        }

        if (sscanf(line, "alpha_Phy1 = %lf", &tmp_val) == 1) {
            alpha[Phy1] = tmp_val;
            biogeo_override_alpha[Phy1] = 1;
        }
        if (sscanf(line, "alpha_Phy2 = %lf", &tmp_val) == 1) {
            alpha[Phy2] = tmp_val;
            biogeo_override_alpha[Phy2] = 1;
        }

        if (sscanf(line, "kmortality_Phy1 = %lf", &tmp_val) == 1) {
            kmortality[Phy1] = tmp_val;
            biogeo_override_kmortality[Phy1] = 1;
        }
        if (sscanf(line, "kmortality_Phy2 = %lf", &tmp_val) == 1) {
            kmortality[Phy2] = tmp_val;
            biogeo_override_kmortality[Phy2] = 1;
        }

        if (sscanf(line, "kox = %lf", &tmp_val) == 1) {
            kox = tmp_val;
            biogeo_override_kox = 1;
        }
        if (sscanf(line, "knit = %lf", &tmp_val) == 1) {
            knit = tmp_val;
            biogeo_override_knit = 1;
        }
        if (sscanf(line, "kcbod_fast = %lf", &tmp_val) == 1) {
            if (isfinite(tmp_val) && tmp_val >= 0.0) kcbod_fast = tmp_val;
        }
        if (sscanf(line, "KO2_cbod = %lf", &tmp_val) == 1) {
            if (isfinite(tmp_val) && tmp_val > 0.0) KO2_cbod = tmp_val;
        }
        if (sscanf(line, "KTOC_fast = %lf", &tmp_val) == 1) {
            if (isfinite(tmp_val) && tmp_val > 0.0) KTOC_fast = tmp_val;
        }
        // Direct wastewater BOD O2 demand
        if (sscanf(line, "k_ww_bod = %lf", &tmp_val) == 1) {
            if (isfinite(tmp_val) && tmp_val >= 0.0) k_ww_bod = tmp_val;
        }
        if (sscanf(line, "KO2_ww_bod = %lf", &tmp_val) == 1) {
            if (isfinite(tmp_val) && tmp_val > 0.0) KO2_ww_bod = tmp_val;
        }
        // Sediment O2 demand
        if (sscanf(line, "SOD_rate = %lf", &tmp_val) == 1) {
            SOD_rate = tmp_val;
        }
        if (sscanf(line, "KO2_SOD = %lf", &tmp_val) == 1) {
            KO2_SOD = tmp_val;
        }
        // Si dissolution fraction
        if (sscanf(line, "ksi_diss = %lf", &tmp_val) == 1) {
            ksi_diss = tmp_val;
        }

        if (sscanf(line, "pCO2atmo = %lf", &tmp_val) == 1) {
            pCO2atmo = tmp_val;
            biogeo_override_pCO2atmo = 1;
        }

        // Gas exchange scaling (applies to both O2 and CO2 piston velocity)
        // Default is 1.0 (no scaling). Values < 1.0 reduce air-water exchange.
        if (sscanf(line, "piston_velocity_scale = %lf", &tmp_val) == 1) {
            if (isfinite(tmp_val) && tmp_val >= 0.0) {
                piston_velocity_scale = tmp_val;
            }
        }

        if (sscanf(line, "KD_Phy = %lf", &tmp_val) == 1) {
            KD_Phy = tmp_val;
        }
        if (sscanf(line, "K_NH4_switch = %lf", &tmp_val) == 1) {
            K_NH4_switch = tmp_val;
        }

        // Mechanistic CO2 limitation for primary production
        if (sscanf(line, "enable_co2_limitation = %d", &enable_co2_limitation) == 1) {}
        if (sscanf(line, "K_CO2_Phy = %lf", &tmp_val) == 1) {
            K_CO2_Phy = tmp_val;
        }

        // Light attenuation and nutrient half-saturation overrides (params.txt)
        if (sscanf(line, "kbg = %lf", &tmp_val) == 1) {
            kbg = tmp_val;
        }
        if (sscanf(line, "kbg_fresh = %lf", &tmp_val) == 1) {
            kbg_fresh = (tmp_val >= 0.0) ? tmp_val : 0.0;
        }
        if (sscanf(line, "S_cdom_threshold = %lf", &tmp_val) == 1) {
            S_cdom_threshold = (tmp_val > 0.0) ? tmp_val : 2.0;
        }
        if (sscanf(line, "S_cdom_width = %lf", &tmp_val) == 1) {
            S_cdom_width = (tmp_val > 0.0) ? tmp_val : 1.0;
        }
        if (sscanf(line, "kbg_upstream = %lf", &tmp_val) == 1) {
            kbg_upstream = (tmp_val >= 0.0) ? tmp_val : 0.0;
        }
        if (sscanf(line, "kbg_transition_km = %lf", &tmp_val) == 1) {
            kbg_transition_km = (tmp_val >= 0.0) ? tmp_val : 0.0;
        }
        if (sscanf(line, "kbg_blend_km = %lf", &tmp_val) == 1) {
            kbg_blend_km = (tmp_val > 0.0) ? tmp_val : 5.0;
        }
        if (sscanf(line, "kspm = %lf", &tmp_val) == 1) {
            kspm = tmp_val;
        }
        if (sscanf(line, "kCDOM = %lf", &tmp_val) == 1) {
            kCDOM = tmp_val;
        }
        if (sscanf(line, "CN_toc = %lf", &tmp_val) == 1) {
            CN_toc = tmp_val;
        }
        if (sscanf(line, "ws_phy = %lf", &tmp_val) == 1) {
            ws_phy = tmp_val;
        }
        if (sscanf(line, "KN_Phy1 = %lf", &tmp_val) == 1) {
            KN[Phy1] = tmp_val;
            biogeo_override_KN[Phy1] = 1;
        }
        if (sscanf(line, "KN_Phy2 = %lf", &tmp_val) == 1) {
            KN[Phy2] = tmp_val;
            biogeo_override_KN[Phy2] = 1;
        }
        if (sscanf(line, "KPO4_Phy1 = %lf", &tmp_val) == 1) {
            KPO4[Phy1] = tmp_val;
            biogeo_override_KPO4[Phy1] = 1;
        }
        if (sscanf(line, "KPO4_Phy2 = %lf", &tmp_val) == 1) {
            KPO4[Phy2] = tmp_val;
            biogeo_override_KPO4[Phy2] = 1;
        }
        if (sscanf(line, "KSi_Phy1 = %lf", &tmp_val) == 1) {
            KSi[Phy1] = tmp_val;
            biogeo_override_KSi = 1;
        }
        if (sscanf(line, "KTOC = %lf", &tmp_val) == 1) {
            if (isfinite(tmp_val) && tmp_val > 0.0) {
                KTOC = tmp_val;
            }
        }
        if (sscanf(line, "kdenit = %lf", &tmp_val) == 1) {
            kdenit = tmp_val;
        }

        // ---------------------------------------------------------------------
        // Phase 5: Carbonate solver failure policy
        // ---------------------------------------------------------------------
        char pol_buf[64];
        if (sscanf(line, "carbonate_solver_failure_policy = %63s", pol_buf) == 1) {
            for (char *c = pol_buf; *c; ++c) {
                *c = (char)toupper((unsigned char)*c);
            }
            if (strcmp(pol_buf, "HOLD_LAST") == 0) {
                carbonate_solver_failure_policy = CARBONATE_FAIL_HOLD_LAST;
            } else if (strcmp(pol_buf, "SKIP_CO2_EXCHANGE") == 0) {
                carbonate_solver_failure_policy = CARBONATE_FAIL_SKIP_CO2_EXCHANGE;
            } else if (strcmp(pol_buf, "FATAL") == 0) {
                carbonate_solver_failure_policy = CARBONATE_FAIL_FATAL;
            } else {
                fprintf(stderr,
                        "Warning: Unknown carbonate_solver_failure_policy '%s'; using HOLD_LAST.\n",
                        pol_buf);
                carbonate_solver_failure_policy = CARBONATE_FAIL_HOLD_LAST;
            }
        }

        if (strstr(line, "riverbed_profile_file =")) {
            if (sscanf(line, "riverbed_profile_file = \"%[^\"]\"", riverbed_profile_file) != 1) {
                fprintf(stderr, "❌ Error: Could not read riverbed_profile_file\n");
            }
        }
        if (strstr(line, "calibration_tidal_targets_file =")) {
            if (sscanf(line, "calibration_tidal_targets_file = \"%[^\"]\"", calibration_tidal_targets_file) != 1) {
                fprintf(stderr, "❌ Error: Could not read calibration_tidal_targets_file\n");
            }
        }

        if (strstr(line, "calibration_parameters_file =")) {
            if (sscanf(line, "calibration_parameters_file = \"%[^\"]\"", calibration_parameters_file) != 1) {
                fprintf(stderr, "❌ Error: Could not read calibration_parameters_file\n");
            }
        }
        if (strstr(line, "calibration_wq_seasonal_targets_file =")) {
            if (sscanf(line, "calibration_wq_seasonal_targets_file = \"%[^\"]\"", calibration_wq_seasonal_targets_file) != 1) {
                fprintf(stderr, "❌ Error: Could not read calibration_wq_seasonal_targets_file\n");
            }
        }

        if (sscanf(line, "num_segments = %d", &num_segments) == 1) {}
        if (sscanf(line, "segment_transition_width = %d", &segment_transition_width) == 1) {}
        if (sscanf(line, "align_segment_break_with_major_tributary = %d", &align_segment_break_with_major_tributary) == 1) {
            align_segment_break_with_major_tributary = align_segment_break_with_major_tributary ? 1 : 0;
        }

        if (sscanf(line, "index_1 = %d", &index_1) == 1) {}
        if (sscanf(line, "B1 = %lf", &B1) == 1) {}
        if (sscanf(line, "LC1 = %lf", &LC1) == 1) {}
        if (sscanf(line, "Chezy1 = %lf", &Chezy1) == 1) {}
        if (sscanf(line, "Rs1 = %lf", &Rs1) == 1) {}

        if (sscanf(line, "index_2 = %d", &index_2) == 1) {}
        if (sscanf(line, "B2 = %lf", &B2) == 1) {}
        if (sscanf(line, "LC2 = %lf", &LC2) == 1) {}
        if (sscanf(line, "Chezy2 = %lf", &Chezy2) == 1) {}
        if (sscanf(line, "Rs2 = %lf", &Rs2) == 1) {}

        if (sscanf(line, "Mero1 = %lf", &Mero1) == 1) {}
        if (sscanf(line, "tau_ero1 = %lf", &tau_ero1) == 1) {}
        if (sscanf(line, "tau_dep1 = %lf", &tau_dep1) == 1) {}

        if (sscanf(line, "Mero2 = %lf", &Mero2) == 1) {}
        if (sscanf(line, "tau_ero2 = %lf", &tau_ero2) == 1) {}
        if (sscanf(line, "tau_dep2 = %lf", &tau_dep2) == 1) {}

        if (sscanf(line, "index_3 = %d", &index_3) == 1) {}
        if (sscanf(line, "Mero3 = %lf", &Mero3) == 1) {}
        if (sscanf(line, "tau_ero3 = %lf", &tau_ero3) == 1) {}
        if (sscanf(line, "tau_dep3 = %lf", &tau_dep3) == 1) {}

        if (sscanf(line, "index_4 = %d", &index_4) == 1) {}
        if (sscanf(line, "Mero4 = %lf", &Mero4) == 1) {}
        if (sscanf(line, "tau_ero4 = %lf", &tau_ero4) == 1) {}
        if (sscanf(line, "tau_dep4 = %lf", &tau_dep4) == 1) {}

        // NOTE: parameter bounds (min/max) are managed by the calibration registry
        // in INPUT/Calibration/calibration_parameters.csv, not by params.txt.

        // Stage-2 P adsorption parameters (used by BIOGEO_P_ADSORPTION_ONLY and full BIOGEO)
        if (sscanf(line, "P_ac = %lf", &P_ac) == 1) {
            if (!isfinite(P_ac) || P_ac <= 0.0) {
                printf("Warning: P_ac invalid; reverting to 35.0 mmol/kg\n");
                P_ac = 35.0;
            }
        }
        if (sscanf(line, "K_ps = %lf", &K_ps) == 1) {
            if (!isfinite(K_ps) || K_ps <= 0.0) {
                printf("Warning: K_ps invalid; reverting to 1.0 mmol/m3\n");
                K_ps = 1.0;
            }
        }
        if (sscanf(line, "k_ads = %lf", &k_ads) == 1) {
            if (!isfinite(k_ads) || k_ads <= 0.0) {
                printf("Warning: k_ads invalid; reverting to 1/3600 s-1\n");
                k_ads = 1.0 / 3600.0;
            }
        }

        // Settling velocity (Winterwerp 2002): ws typically 0.1-2.0 mm/s for estuarine flocs
        if (sscanf(line, "settling_velocity = %lf", &ws) == 1) {
            if (!isfinite(ws) || ws <= 0.0 || ws > 0.01) {
                printf("Warning: settling_velocity invalid (must be 0 < ws <= 0.01 m/s); reverting to 5e-4 m/s\n");
                ws = 5.0e-4;  // Default 0.5 mm/s
            }
        }

        // C_VDB is a SCALING FACTOR for the Van der Burgh K predictor, not K itself
        // Range: 0.1-2.0 to achieve literature K values of 0.2-0.8 (Savenije 2012 §4.4)
        if (sscanf(line, "C_VDB = %lf", &C_VDB) == 1) {
            if (!isfinite(C_VDB) || C_VDB <= 0.0) {
                printf("Warning: C_VDB invalid; reverting to 0.5\n");
                C_VDB = 0.5;
            } else if (C_VDB < 0.1 || C_VDB > 2.0) {
                printf("Warning: C_VDB %.3f outside typical range [0.1-2.0], proceeding with caution\n", C_VDB);
            }
        }

        if (sscanf(line, "D0_CORRECTION = %lf", &D0_CORRECTION) == 1) {
            if (!isfinite(D0_CORRECTION) || D0_CORRECTION <= 0.0) {
                printf("Warning: D0_CORRECTION invalid; reverting to 1.0\n");
                D0_CORRECTION = 1.0;
            } else if (D0_CORRECTION < 0.3 || D0_CORRECTION > 3.0) {
                printf("Warning: D0_CORRECTION %.3f outside typical range [0.3-3.0], proceeding with caution\n", D0_CORRECTION);
            }
        }

        if (sscanf(line, "K_DISCHARGE_SENSITIVITY = %lf", &K_DISCHARGE_SENSITIVITY) == 1) {
            if (!isfinite(K_DISCHARGE_SENSITIVITY) || K_DISCHARGE_SENSITIVITY < 0.0) {
                printf("Warning: K_DISCHARGE_SENSITIVITY invalid; reverting to 0.0\n");
                K_DISCHARGE_SENSITIVITY = 0.0;
            } else if (K_DISCHARGE_SENSITIVITY > 2.0) {
                printf("Warning: K_DISCHARGE_SENSITIVITY %.3f > 2.0, proceeding with caution\n", K_DISCHARGE_SENSITIVITY);
            }
        }

        if (sscanf(line, "dispersion_abs_min = %lf", &dispersion_abs_min) == 1) {}
        if (sscanf(line, "dispersion_abs_max = %lf", &dispersion_abs_max) == 1) {}

        // Discharge-driven turbidity parameters
        if (sscanf(line, "k_turb_Q = %lf", &k_turb_Q) == 1) {
            if (!isfinite(k_turb_Q) || k_turb_Q < 0.0) {
                k_turb_Q = 0.0;
            }
        }
        if (sscanf(line, "Q_turb_ref = %lf", &Q_turb_ref) == 1) {
            if (!isfinite(Q_turb_ref) || Q_turb_ref < 0.0) {
                Q_turb_ref = 0.0;  // 0 = auto-compute from forcing data
            }
        }
        if (sscanf(line, "alpha_turb = %lf", &alpha_turb) == 1) {
            if (!isfinite(alpha_turb) || alpha_turb < 0.0 || alpha_turb > 5.0) {
                alpha_turb = 1.5;
            }
        }
        // Depth-dependent resuspension turbidity
        if (sscanf(line, "k_resus = %lf", &k_resus) == 1) {
            if (!isfinite(k_resus) || k_resus < 0.0) {
                k_resus = 0.0;
            }
        }
    }

    fclose(file);

    // ------------------------------------------------------------------
    // Strict scientific mode for reproducibility
    // ------------------------------------------------------------------
    // No environment-based physics/runtime overrides or hidden runtime caps.
    // Runtime clocks and physics must be fully versioned in INPUT/params.txt.
    if (MAXT <= WARMUP) {
        fprintf(stderr,
                "❌ Error: MAXT (%ld days) must exceed WARMUP (%ld days) in INPUT/params.txt.\n",
                MAXT, WARMUP);
        exit(EXIT_FAILURE);
    }

    MAXT *= (24 * 60 * 60);
    WARMUP *= (24 * 60 * 60);

    long main_window_seconds = MAXT - WARMUP;
    if (main_window_seconds <= 0) {
        fprintf(stderr, "❌ Error: MAXT must exceed WARMUP to define a diagnostic window.\n");
        exit(EXIT_FAILURE);
    }

    simulation_window_seconds = main_window_seconds;

    PRINTF_DATA("   MAXT = %ld, WARMUP = %ld, DELTI = %d, TS = %d\n", MAXT, WARMUP, DELTI, TS);
    PRINTF_DATA("   DELXI = %d, EL = %d, AMPL = %.2f, num_segments = %d, trans_width = %d\n", DELXI, EL, AMPL, num_segments, segment_transition_width);
    PRINTF_DATA("   Segment 1: index_1 = %d, B1 = %.2f, LC1 = %.2f, Chezy1 = %.2f, Rs1 = %.2f\n",
                index_1, B1, LC1, Chezy1, Rs1);
    PRINTF_DATA("   Segment 2: index_2 = %d, B2 = %.2f, LC2 = %.2f, Chezy2 = %.2f, Rs2 = %.2f\n",
                index_2, B2, LC2, Chezy2, Rs2);
    PRINTF_DATA("   Mero1 = %.6f, tau_ero1 = %.2f, tau_dep1 = %.2f\n", Mero1, tau_ero1, tau_dep1);
    PRINTF_DATA("   Mero2 = %.6f, tau_ero2 = %.2f, tau_dep2 = %.2f\n", Mero2, tau_ero2, tau_dep2);

    M = EL / DELXI;
    M = (M % 2 == 0) ? M : M + 1;

    if (M >= MAXM - 1) {
        printf("WARNING: Computed M=%d exceeds safe array bounds. Limiting to %d\n", M, MAXM - 2);
        M = MAXM - 2;
    }

    M1 = M - 1;
    M2 = M - 2;
    M3 = M - 3;

    if (calibration_mode == 0) {
        printf("✅ Computed: M = %d, M1 = %d, M2 = %d, M3 = %d\n", M, M1, M2, M3);
    }

    if (M <= 0) {
        printf("❌ ERROR: Computed M is non-positive! Exiting...\n");
        exit(EXIT_FAILURE);
    }
}

// -----------------------------------------------------------------------------
// Riverbed profile and configuration loaders
// -----------------------------------------------------------------------------

void readRiverbedProfile(double *riverbed_depth, const char *filename) {
    double raw_depths[MAXM + 1];
    memset(raw_depths, 0, sizeof(raw_depths));
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "❌ Error: Cannot open riverbed profile file %s\n", filename);
        exit(1);
    }

    char line[256];
    int cell;
    double depth;
    int max_read_cell = 0;

    if (fgets(line, sizeof(line), fp) && line[0] == '#') {
        printf("✅ Skipping header line: %s", line);
    } else {
        rewind(fp);
    }

    int lineNum = 0;
    double max_depth = 0.0;
    double min_depth = 1000.0;
    double sum_depth = 0.0;
    int valid_points = 0;

    while (fgets(line, sizeof(line), fp) != NULL) {
        lineNum++;
        if (sscanf(line, "%d,%lf", &cell, &depth) == 2) {
            if (cell >= 0 && cell <= M) {
                if (depth < 0.0 || isnan(depth)) {
                    fprintf(stderr, "⚠️ Warning: Invalid depth %.2f for cell %d. Setting to 0.\n",
                            depth, cell);
                    depth = 0.0;
                } else {
                    max_depth = (depth > max_depth) ? depth : max_depth;
                    min_depth = (depth < min_depth) ? depth : min_depth;
                    sum_depth += depth;
                    valid_points++;
                }

                raw_depths[cell] = depth;
                max_read_cell = (cell > max_read_cell) ? cell : max_read_cell;
            } else {
                fprintf(stderr, "⚠️ Warning: Cell index %d out of range (0-%d) on line %d. Skipping.\n",
                        cell, M, lineNum);
            }
        } else {
            fprintf(stderr, "⚠️ Warning: Could not parse line %d: '%s'. Skipping.\n", lineNum, line);
        }
    }

    fclose(fp);

    if (max_read_cell < M) {
        PRINTF_DATA("⚠️ Warning: Riverbed profile data incomplete. Filling missing values with last known depth.\n");
        for (int i = max_read_cell + 1; i <= M; i++) {
            raw_depths[i] = raw_depths[max_read_cell];
        }
    }

    double avg_depth = valid_points > 0 ? sum_depth / valid_points : 0.0;
    double max_allowed_jump = avg_depth * 0.3;

    PRINTF_DATA("🔄 Applying smoothing to riverbed profile to reduce numerical instabilities...\n");
    int window_size = 3;
    int half_window = window_size / 2;

    for (int i = 0; i <= M; i++) {
        if (i < half_window || i > M - half_window) {
            riverbed_depth[i] = raw_depths[i];
        } else {
            double sum = 0.0;
            for (int j = -half_window; j <= half_window; j++) {
                sum += raw_depths[i + j];
            }
            riverbed_depth[i] = sum / window_size;
        }
    }

    int abrupt_changes = 0;
    for (int i = 1; i <= M; i++) {
        if (i > 1 && fabs(riverbed_depth[i] - riverbed_depth[i - 1]) > max_allowed_jump) {
            abrupt_changes++;
            printf("⚠️ Smoothing large depth jump between cells %d and %d (%.2f → %.2f)\n",
                   i - 1, i, riverbed_depth[i - 1], riverbed_depth[i]);

            double depth_diff = riverbed_depth[i] - riverbed_depth[i - 1];
            int smooth_window = 5;

            if (i + smooth_window <= M) {
                double original_value = riverbed_depth[i];
                riverbed_depth[i] = riverbed_depth[i - 1] + 0.25 * depth_diff;
                for (int j = 1; j < smooth_window; j++) {
                    double blend = (double)(j + 1) / smooth_window;
                    riverbed_depth[i + j] = riverbed_depth[i - 1] + blend * depth_diff;
                }
                printf("  → Smoothed transition from %.2f to %.2f over %d cells\n",
                       riverbed_depth[i - 1], original_value, smooth_window);
            }
        }
    }

    if (calibration_mode == 0) {
        printf("✅ Riverbed profile processed: %d points, depth range: %.2f-%.2f m, avg: %.2f m\n",
               valid_points, min_depth, max_depth, avg_depth);
    }
    if (abrupt_changes > 0) {
        printf("✅ Smoothed %d abrupt changes in riverbed profile\n", abrupt_changes);
    }
}

void readConfigFile(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("❌ Error opening config file");
        exit(1);
    }

    char line[256];
    int foundCount = 0;
    int missingCount = 0;
    int boundaryType = -1;

    while (fgets(line, sizeof(line), fp)) {
        trim_line_inplace(line);

        if (line[0] == '#' || line[0] == '\0') continue;

        if (strstr(line, "type=")) {
            if (strstr(line, "UpperBoundary")) {
                boundaryType = 0;
            } else if (strstr(line, "LowerBoundary")) {
                boundaryType = 1;
            } else if (strstr(line, "Forcing")) {
                boundaryType = 2;
            } else if (strstr(line, "Tributary")) {
                boundaryType = -1;
                continue;
            }
        }

        if (boundaryType < 0) continue;

        int found = 0;

        if (boundaryType == 2) {
            for (int f = 0; f < FORCING_COUNT; f++) {
                char key[64];
                sprintf(key, "%sFile", forcingNames[f]);
                if (strstr(line, key)) {
                    forcingData[f].filePath = strdup(strchr(line, '=') + 1);
                    forcingData[f].required = 1;
                    found = 1;
                    foundCount++;
                    break;
                }
            }
            continue;
        }

        for (int i = 0; i < CHEM_COUNT; i++) {
            char key[64];
            sprintf(key, "%s=", variableNames[i]);

            if (strstr(line, key) == line) {
                if (boundaryType == 0) {
                    upstreamBC[i].filePath = strdup(strchr(line, '=') + 1);
                    upstreamBC[i].required = isDataInputRequired(i);
                } else {
                    downstreamBC[i].filePath = strdup(strchr(line, '=') + 1);
                    downstreamBC[i].required = isDataInputRequired(i);
                }
                found = 1;
                foundCount++;
                break;
            }
        }

        if (!found && boundaryType == 2) {
            if (strstr(line, "DischargeFile=")) {
                forcingData[FORCING_DISCHARGE].filePath = strdup(strchr(line, '=') + 1);
                foundCount++;
                found = 1;
            } else if (strstr(line, "ElevationFile=")) {
                forcingData[FORCING_ELEVATION].filePath = strdup(strchr(line, '=') + 1);
                foundCount++;
                found = 1;
            } else if (strstr(line, "LightFile=")) {
                forcingData[FORCING_LIGHT].filePath = strdup(strchr(line, '=') + 1);
                foundCount++;
                found = 1;
            } else if (strstr(line, "TemperatureFile=")) {
                forcingData[FORCING_TEMPERATURE].filePath = strdup(strchr(line, '=') + 1);
                foundCount++;
                found = 1;
            } else if (strstr(line, "WindFile=")) {
                forcingData[FORCING_WIND].filePath = strdup(strchr(line, '=') + 1);
                foundCount++;
                found = 1;
            }
        }

        if (!found && strstr(line, "riverbed_profile_file=")) {
            strcpy(riverbed_profile_file, strchr(line, '=') + 1);
            foundCount++;
        }
    }

    fclose(fp);

    for (int f = 0; f < FORCING_COUNT; f++) {
        if (!isValidPath(forcingData[f].filePath)) {
            printf("❌ Error: Required forcing data '%s' is missing!\n", forcingNames[f]);
            missingCount++;
        }
    }

    checkFilePath(riverbed_profile_file, "riverbed_profile_file");

    if (missingCount > 0) {
        fprintf(stderr, "❌ Error: %d required input files are missing. Cannot continue.\n", missingCount);
        exit(EXIT_FAILURE);
    }
}

static int validateChemicalBoundaryData(void);
static int validateForcingData(void);

void readBoundaryData(const char *filename) {
    PRINTF_DATA("📂 Reading boundary conditions from: %s\n", filename);

    readConfigFile(filename);

    int size = MAX_TIMESERIES;

    for (int i = 0; i < CHEM_COUNT; i++) {
        upstreamBC[i].time = (double *)malloc(size * sizeof(double));
        upstreamBC[i].data = (double *)malloc(size * sizeof(double));
        if (upstreamBC[i].time == NULL || upstreamBC[i].data == NULL) {
            fprintf(stderr, "❌ Error allocating memory for upstream BC: %s\n", variableNames[i]);
            exit(EXIT_FAILURE);
        }
    }

    for (int i = 0; i < CHEM_COUNT; i++) {
        downstreamBC[i].time = (double *)malloc(size * sizeof(double));
        downstreamBC[i].data = (double *)malloc(size * sizeof(double));
        if (downstreamBC[i].time == NULL || downstreamBC[i].data == NULL) {
            fprintf(stderr, "❌ Error allocating memory for downstream BC: %s\n", variableNames[i]);
            exit(EXIT_FAILURE);
        }
    }

    for (int f = 0; f < FORCING_COUNT; f++) {
        forcingData[f].time = (double *)malloc(size * sizeof(double));
        forcingData[f].data = (double *)malloc(size * sizeof(double));
        if (forcingData[f].time == NULL || forcingData[f].data == NULL) {
            fprintf(stderr, "❌ Error allocating memory for forcing data: %s\n", forcingNames[f]);
            exit(EXIT_FAILURE);
        }
    }

    for (int i = 0; i < CHEM_COUNT; i++) {
        if (isValidPath(upstreamBC[i].filePath)) {
            upstreamBC[i].dataSize = readFile(MAX_TIMESERIES, upstreamBC[i].time, upstreamBC[i].data, upstreamBC[i].filePath);
            if (upstreamBC[i].dataSize <= 0) {
                printf("⚠️ Warning: No data read for upstream boundary '%s'\n", variableNames[i]);
            }
        } else if (upstreamBC[i].required && i != Sal) {
            printf("❌ Error: Required upstream boundary '%s' file not found or invalid\n", variableNames[i]);
            upstreamBC[i].dataSize = 0;
        }
    }

    for (int i = 0; i < CHEM_COUNT; i++) {
        if (isValidPath(downstreamBC[i].filePath)) {
            downstreamBC[i].dataSize = readFile(MAX_TIMESERIES, downstreamBC[i].time, downstreamBC[i].data, downstreamBC[i].filePath);
            if (downstreamBC[i].dataSize <= 0) {
                printf("⚠️ Warning: No data read for downstream boundary '%s'\n", variableNames[i]);
            }
        } else if (downstreamBC[i].required) {
            printf("❌ Error: Required downstream boundary '%s' file not found or invalid\n", variableNames[i]);
            downstreamBC[i].dataSize = 0;
        }
    }

    for (int f = 0; f < FORCING_COUNT; f++) {
        if (isValidPath(forcingData[f].filePath)) {
            forcingData[f].dataSize = readFile(MAX_TIMESERIES, forcingData[f].time,
                                             forcingData[f].data,
                                             forcingData[f].filePath);
            if (forcingData[f].dataSize <= 0) {
                fprintf(stderr, "❌ Error: No data read for forcing: %s\n",
                        forcingNames[f]);
                exit(EXIT_FAILURE);
            }
        } else if (forcingData[f].required) {
            printf("❌ Error: Required forcing data '%s' file not found or invalid\n", forcingNames[f]);
        }
    }

    setupChemicalBoundaryDataMap();

    if (!validateChemicalBoundaryData() || !validateForcingData()) {
        printf("⚠️  Error: Update config_input.txt for missing files!\n");
        exit(EXIT_FAILURE);
    }
}

void readGlobalConfigSettings(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("⚠️ Warning: Cannot open config file %s for global settings. Using defaults.\n", filename);
        tributaryEnabled = 0;
        return;
    }

    char line[256];

    while (fgets(line, sizeof(line), fp)) {
        trim_line_inplace(line);

        if (line[0] == '#' || line[0] == '\0') continue;

        if (strstr(line, "tributaryEnabled=") == line) {
            int value = 0;
            if (sscanf(line, "tributaryEnabled=%d", &value) == 1) {
                tributaryEnabled = (value != 0) ? 1 : 0;
                PRINTF_DATA("✅ Tributary system %s (tributaryEnabled=%d)\n",
                            tributaryEnabled ? "ENABLED" : "DISABLED", tributaryEnabled);
            }
        }

        if (strstr(line, "name=") || strstr(line, "type=")) {
            break;
        }
    }

    fclose(fp);
}

void readTributaryData(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("⚠️ Warning: Cannot open config file %s for tributary data. No tributaries will be used.\n", filename);
        numTributaries = 0;
        return;
    }

    PRINTF_DATA("📂 Reading tributary configuration from: %s\n", filename);

    char line[256];

    numTributaries = 0;
    while (fgets(line, sizeof(line), fp)) {
        trim_line_inplace(line);

        if (line[0] == '#' || line[0] == '\0') continue;

        if (strstr(line, "type=") && strstr(line, "Tributary")) {
            numTributaries++;
        }
    }

    if (numTributaries <= 0) {
        fclose(fp);
        printf("ℹ️ No tributaries configured\n");
        return;
    }

    tributaries = (Tributary *)malloc(numTributaries * sizeof(Tributary));
    if (!tributaries) {
        perror("❌ Memory allocation failed for tributaries");
        exit(1);
    }

    for (int i = 0; i < numTributaries; i++) {
        memset(&tributaries[i], 0, sizeof(Tributary));
        sprintf(tributaries[i].name, "Trib%d", i + 1);
        for (int c = 0; c < CHEM_COUNT; c++) {
            tributaries[i].chemicalData[c].filePath[0] = '\0';
            tributaries[i].chemicalData[c].timeArray = NULL;
            tributaries[i].chemicalData[c].dataArray = NULL;
        }
    }

    rewind(fp);
    int currentTrib = -1;

    while (fgets(line, sizeof(line), fp)) {
        trim_line_inplace(line);

        if (line[0] == '#' || line[0] == '\0') continue;

        if (strstr(line, "name=") == line) {
            char name[50] = {0};
            sscanf(line, "name=%s", name);

            char nextLine[256];
            long pos = ftell(fp);
            if (fgets(nextLine, sizeof(nextLine), fp)) {
                if (strstr(nextLine, "type=") && strstr(nextLine, "Tributary")) {
                    currentTrib++;
                    if (currentTrib < numTributaries) {
                        strncpy(tributaries[currentTrib].name, name, sizeof(tributaries[currentTrib].name) - 1);
                    }
                }
            }
            fseek(fp, pos, SEEK_SET);
        }

        if (currentTrib >= 0 && currentTrib < numTributaries) {
            if (strstr(line, "cellIndex=") == line) {
                sscanf(line, "cellIndex=%d", &tributaries[currentTrib].cellIndex);
            } else if (strstr(line, "is_urban=") == line) {
                sscanf(line, "is_urban=%d", &tributaries[currentTrib].is_urban);
            } else if (strstr(line, "discharge=") == line) {
                const char *path = strchr(line, '=') + 1;
                strncpy(tributaries[currentTrib].dischargeFile, path, sizeof(tributaries[currentTrib].dischargeFile) - 1);
            } else {
                for (int c = 0; c < CHEM_COUNT; c++) {
                    char varPattern[64];
                    sprintf(varPattern, "%s=", variableNames[c]);

                    if (strstr(line, varPattern) == line) {
                        const char *path = strchr(line, '=') + 1;
                        strncpy(tributaries[currentTrib].chemicalData[c].filePath, path,
                                sizeof(tributaries[currentTrib].chemicalData[c].filePath) - 1);
                        break;
                    }
                }
            }
        }
    }

    fclose(fp);

    // Strict scientific mode: tributary placement is controlled only by INPUT/config_input.txt.

    for (int j = 0; j < numTributaries; j++) {
        PRINTF_DATA("📊 Tributary %d: %s at cell index %d\n", j + 1,
                    tributaries[j].name, tributaries[j].cellIndex);

        if (tributaries[j].cellIndex < 0 || tributaries[j].cellIndex > M) {
            PRINTF_ERROR("⚠️ Warning: Tributary %s has invalid cell index %d (valid range: 0-%d)\n",
                         tributaries[j].name, tributaries[j].cellIndex, M);
        }

        if (!isValidPath(tributaries[j].dischargeFile)) {
            printf("❌ Error: Required parameter 'discharge' missing for tributary %s\n",
                   tributaries[j].name);
            exit(EXIT_FAILURE);
        }
    }

    int size = MAX_TIMESERIES;

    for (int j = 0; j < numTributaries; j++) {
        tributaries[j].dischargeTime = (double *)malloc(size * sizeof(double));
        tributaries[j].discharge = (double *)malloc(size * sizeof(double));

        if (!tributaries[j].dischargeTime || !tributaries[j].discharge) {
            fprintf(stderr, "❌ Error allocating memory for tributary NET discharge\n");
            exit(EXIT_FAILURE);
        }

        int actualDataSize = readFile(size, tributaries[j].dischargeTime, tributaries[j].discharge,
                                     tributaries[j].dischargeFile);

        tributaries[j].dischargeDataSize = actualDataSize;

        if (debug_level > 0) {
            PRINTF_DATA("📊 Tributary %s: Loaded %d discharge data points\n",
                        tributaries[j].name, actualDataSize);
        }

        for (int c = 0; c < CHEM_COUNT; c++) {
            if (c == pCO2 || c == PH || c == CO2) continue;

            if (isValidPath(tributaries[j].chemicalData[c].filePath)) {
                tributaries[j].chemicalData[c].timeArray = (double *)malloc(size * sizeof(double));
                tributaries[j].chemicalData[c].dataArray = (double *)malloc(size * sizeof(double));

                if (!tributaries[j].chemicalData[c].timeArray ||
                    !tributaries[j].chemicalData[c].dataArray) {
                    fprintf(stderr, "❌ Error allocating memory for tributary %s variable %s\n",
                            tributaries[j].name, variableNames[c]);
                    exit(EXIT_FAILURE);
                }

                int nread = readFile(size, tributaries[j].chemicalData[c].timeArray,
                                     tributaries[j].chemicalData[c].dataArray,
                                     tributaries[j].chemicalData[c].filePath);
                tributaries[j].chemicalData[c].dataSize = nread;

                // Transparency guard: flag near-constant tributary chemistry series.
                // Constant forcings may be intentional for data-sparse cases, but
                // should be visible to users/reviewers during calibration.
                if (nread > 30) {
                    double vmin = tributaries[j].chemicalData[c].dataArray[0];
                    double vmax = vmin;
                    for (int k = 1; k < nread; ++k) {
                        double vv = tributaries[j].chemicalData[c].dataArray[k];
                        if (vv < vmin) vmin = vv;
                        if (vv > vmax) vmax = vv;
                    }
                    double span = vmax - vmin;
                    if (isfinite(span) && span <= 1e-9) {
                        printf("⚠️  Tributary %s %s appears constant (n=%d, value=%.6g). Verify data basis.\n",
                               tributaries[j].name, variableNames[c], nread, vmin);
                    }
                }
            }
        }
    }
}

static int validateChemicalBoundaryData() {
    int valid = 1;
    int upstream_count = 0;
    int downstream_count = 0;

    for (int var = 0; var < NUM_BC_VARS; var++) {
        if (var == pCO2 || var == PH || var == CO2) continue;

        if (isDataInputRequired(var)) {
            if (var == Sal) {
                if (isValidPath(downstreamBC[var].filePath)) {
                    downstream_count++;
                } else {
                    printf("❌ Error: Missing downstream boundary data for %s\n", variableNames[var]);
                    valid = 0;
                }
                continue;
            }

            if (isValidPath(upstreamBC[var].filePath)) {
                upstream_count++;
            } else {
                printf("❌ Error: Missing upstream boundary data for %s\n", variableNames[var]);
                valid = 0;
            }

            if (isValidPath(downstreamBC[var].filePath)) {
                downstream_count++;
            } else {
                printf("❌ Error: Missing downstream boundary data for %s\n", variableNames[var]);
                valid = 0;
            }
        }
    }

    if (upstream_count == 0 || downstream_count == 0) {
        printf("❌ Error: No chemical boundary data found for %s boundary!\n",
               upstream_count == 0 ? "upstream" : "downstream");
        valid = 0;
    }

    return valid;
}

static int validateForcingData() {
    int valid = 1;

    for (int f = 0; f < FORCING_COUNT; f++) {
        if (forcingData[f].required) {
            if (!forcingData[f].time || !forcingData[f].data) {
                printf("❌ Error: Missing required forcing data arrays for %s\n",
                       forcingNames[f]);
                valid = 0;
            }
            if (forcingData[f].dataSize <= 0) {
                printf("❌ Error: Invalid data size (%d) for forcing data %s\n",
                       forcingData[f].dataSize, forcingNames[f]);
                valid = 0;
            }
        }
    }

    return valid;
}
