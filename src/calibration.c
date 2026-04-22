/**
 * @file calibration.c
 * @brief Stage-based calibration using unified long-form seasonal targets
 *
 * A single seasonal target file can contain rows for stages 1..4 via a per-row
 * stage column; only rows matching the active calibration stage are scored.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>

#ifdef _WIN32
#include <direct.h>
#define strcasecmp _stricmp
#endif

#include <nlopt.h>

#include "define.h"
#include "variables.h"
#include "file.h"
#include "diagnostics.h"
#include "utilities.h"
#include "init.h"
#include "transport.h"
#include "biogeo.h"
#include "calibration_helpers.h"

// Stage sizes
#define NUM_PARAMS_STAGE1 9  // Hydro + Transport
#define NUM_PARAMS_STAGE2 15 // Sediment (4 segments) + P-Adsorption
#define NUM_PARAMS_STAGE3 20 // Eutrophication + light + nutrient half-saturations + O2/SOD + quadratic closure + fast CBOD + CDOM + CN_toc
#define NUM_PARAMS_STAGE4 1  // Carbonate (Piston velocity)

// Stage 1: Hydrodynamics (Chezy, Rs) + Transport (Dispersion)
enum {
    IDX1_CHEZY1 = 0,
    IDX1_CHEZY2,
    IDX1_LC1,
    IDX1_LC2,
    IDX1_RS1,
    IDX1_RS2,
    IDX1_C_VDB,
    IDX1_D0_CORR,
    IDX1_AMPL_CORR,
};

// Stage 2: Sediment (Mero, tau) + Phosphorus (Langmuir)
enum {
    IDX2_MERO1 = 0,
    IDX2_TAU_ERO1,
    IDX2_TAU_DEP1,
    IDX2_MERO2,
    IDX2_TAU_ERO2,
    IDX2_TAU_DEP2,
    IDX2_MERO3,
    IDX2_TAU_ERO3,
    IDX2_TAU_DEP3,
    IDX2_MERO4,
    IDX2_TAU_ERO4,
    IDX2_TAU_DEP4,
    IDX2_P_AC,
    IDX2_K_PS,
    IDX2_K_ADS,
};

// Stage 3: Eutrophication (Rate multipliers)
enum {
    IDX3_SCALE_KOX = 0,
    IDX3_SCALE_KNIT,
    IDX3_SCALE_KDENIT,
    IDX3_SCALE_PB_PHY,     // Unified max production multiplier
    IDX3_SCALE_KMORT_PHY,  // Unified mortality multiplier
    IDX3_SCALE_KD_PHY,
    IDX3_SCALE_K_NH4_SWITCH,
    IDX3_SCALE_KBG,      // Background light attenuation multiplier
    IDX3_SCALE_KSPM,     // SPM light attenuation multiplier
    IDX3_SCALE_WS_PHY,   // Phytoplankton settling velocity multiplier
    IDX3_SCALE_ALPHA,    // Unified light efficiency multiplier
    IDX3_SCALE_KN,       // Nitrate half-saturation multiplier
    IDX3_SCALE_KPO4,     // Phosphate half-saturation multiplier
    IDX3_SCALE_KSI,      // Silica half-saturation multiplier
    IDX3_SCALE_KTOC,     // TOC half-saturation multiplier
    IDX3_SOD_RATE,       // Direct SOD calibration [mmol O2/m2/s]
    IDX3_KMORT2_PHY1,     // Direct quadratic mortality closure [s^-1/(mmolC/m^3)]
    IDX3_KCBOD_FAST,      // Direct fast labile-C degradation rate [s^-1]
    IDX3_KCDOM,           // CDOM light absorption [m^-1/(mgC/L)]
    IDX3_CN_TOC,          // Effective C:N for bulk TOC mineralization [molC/molN]
};

// Stage 4: Carbonate
enum {
    IDX4_SCALE_PISTON_VEL = 0,
};

typedef enum {
    PARAM_TRANSFORM_LINEAR = 0,
    PARAM_TRANSFORM_LOG10 = 1,
} ParamTransform;

// Strict CSV parsing helpers (no quoting support; canonical files must be clean)

static void trim_trailing_whitespace(char *str)
{
    if (!str) return;
    size_t n = strlen(str);
    while (n > 0 && (str[n - 1] == '\n' || str[n - 1] == '\r' || isspace((unsigned char)str[n - 1]))) {
        str[n - 1] = '\0';
        n--;
    }
}

static void strip_utf8_bom(char *s)
{
    if (!s) return;
    const unsigned char *u = (const unsigned char *)s;
    if (u[0] == 0xEF && u[1] == 0xBB && u[2] == 0xBF) {
        size_t len = strlen(s);
        if (len >= 3) {
            memmove(s, s + 3, len - 2);
        }
    }
}

static bool try_parse_double_strict(const char *text, double *value_out)
{
    if (!text || !value_out) return false;
    while (isspace((unsigned char)*text)) text++;

    char *end = NULL;
    errno = 0;
    double val_double = strtod(text, &end);
    if (errno != 0 || end == text) return false;

    while (end && *end != '\0' && isspace((unsigned char)*end)) end++;
    if (end && *end != '\0') return false;

    if (!isfinite(val_double)) return false;
    *value_out = val_double;
    return true;
}

static bool try_parse_int_strict(const char *text, int *value_out)
{
    if (!text || !value_out) return false;
    while (isspace((unsigned char)*text)) text++;

    char *end = NULL;
    errno = 0;
    long val_int = strtol(text, &end, 10);
    if (errno != 0 || end == text) return false;

    while (end && *end != '\0' && isspace((unsigned char)*end)) end++;
    if (end && *end != '\0') return false;

    if (val_int < INT_MIN || val_int > INT_MAX) return false;
    *value_out = (int)val_int;
    return true;
}

static bool try_parse_enabled_flag(const char *text, int *value_out)
{
    if (!text || !value_out) return false;

    int ival = 0;
    if (try_parse_int_strict(text, &ival)) {
        *value_out = (ival != 0) ? 1 : 0;
        return true;
    }

    double dval = NAN;
    if (try_parse_double_strict(text, &dval)) {
        *value_out = (fabs(dval) > 0.0) ? 1 : 0;
        return true;
    }

    // Accept common textual booleans for robustness.
    if (strcasecmp(text, "true") == 0 || strcasecmp(text, "yes") == 0 || strcasecmp(text, "on") == 0) {
        *value_out = 1;
        return true;
    }
    if (strcasecmp(text, "false") == 0 || strcasecmp(text, "no") == 0 || strcasecmp(text, "off") == 0) {
        *value_out = 0;
        return true;
    }

    return false;
}

static bool split_csv_row_n(char *line, char **tokens_out, int ncols)
{
    // Expects exactly ncols columns and rejects extra columns.
    int token_count = 0;
    char *token = strtok(line, ",");
    while (token && token_count < ncols) {
        tokens_out[token_count++] = token;
        token = strtok(NULL, ",");
    }

    if (token != NULL) return false; // >ncols columns
    return token_count == ncols;
}

// Target structures

typedef struct {
    double location_km;
    double range_m;
} TidalTarget;

typedef struct {
    TidalTarget *points;
    size_t count;
} TidalTargetSet;

typedef enum {
    WQ_VAR_CHEM = 0,
    WQ_VAR_PHY_TOTAL = 1,
} WQVarKind;

static const double CHLA_MOLAR_MASS_C_MG_PER_MMOL = 12.01;

typedef struct {
    // Timestamp convention: seconds since post-warmup (analysis) start.
    // I.e., targets are defined on the t_rel = (t_sim - WARMUP) time base.
    long target_step_seconds;
    int cell_i0;
    int cell_i1;
    double alpha;
    int stage; // 1..4
    WQVarKind kind;
    Chem chem; // valid if kind==WQ_VAR_CHEM
    bool is_chla_ugl; // if true and kind==WQ_VAR_PHY_TOTAL, model value uses Phy->Chla observation operator
    double target_value;
    double scale;
    double weight;

    bool recorded;
    double model_value;
} SeasonalWQTarget;

typedef struct {
    SeasonalWQTarget *points;
    size_t count;
} SeasonalWQTargetSet;

static TidalTargetSet g_tidal = { NULL, 0 };
static SeasonalWQTargetSet g_wq_seasonal = { NULL, 0 };

static const char *chem_short_name(Chem c)
{
    switch (c) {
        case Phy1: return "Phy1";
        case Phy2: return "Phy2";
        case Si: return "Si";
        case NO3: return "NO3";
        case NH4: return "NH4";
        case PO4: return "PO4";
        case PIP: return "PIP";
        case O2: return "O2";
        case TOC: return "TOC";
        case Sal: return "Sal";
        case SPM: return "SPM";
        case DIC: return "DIC";
        case AT: return "AT";
        case pCO2: return "pCO2";
        case PH: return "pH";
        case CO2: return "CO2";
        default: return "?";
    }
}

static void print_stage_variable_counts_inline(const CGEM_WQStageSummary *s)
{
    if (!s) return;
    int any = 0;
    if (s->phy_total > 0) {
        printf(" Phy(total)=%zu", s->phy_total);
        any = 1;
    }
    for (int c = 0; c < CHEM_COUNT; ++c) {
        if (s->by_chem[c] > 0) {
            printf(" %s=%zu", chem_short_name((Chem)c), s->by_chem[c]);
            any = 1;
        }
    }
    if (!any) printf(" (none)");
}

static void log_loaded_seasonal_targets_summary(int active_stage)
{
    if (g_wq_seasonal.count == 0) return;

    CGEM_WQTargetSummaryInput *in = (CGEM_WQTargetSummaryInput *)malloc(g_wq_seasonal.count * sizeof(CGEM_WQTargetSummaryInput));
    if (!in) return;

    for (size_t i = 0; i < g_wq_seasonal.count; ++i) {
        const SeasonalWQTarget *p = &g_wq_seasonal.points[i];
        in[i].target_step_seconds = p->target_step_seconds;
        in[i].stage = p->stage;
        in[i].kind = (p->kind == WQ_VAR_PHY_TOTAL) ? CGEM_WQ_TARGET_KIND_PHY_TOTAL : CGEM_WQ_TARGET_KIND_CHEM;
        in[i].chem = p->chem;
    }

    CGEM_WQTargetSummary sum;
    if (cgem_compute_wq_target_summary(in, g_wq_seasonal.count, &sum) != 0) {
        free(in);
        return;
    }
    free(in);

    if (debug_level >= DEBUG_LEVEL_ESSENTIAL) {
        printf("\n📌 Seasonal targets loaded: total=%zu (active_stage=%d)\n", sum.total_targets, active_stage);
        for (int st = 1; st <= 4; ++st) {
            const CGEM_WQStageSummary *s = &sum.stage[st];
            if (s->total == 0) continue;
             printf("  Stage %d: n=%zu, step_seconds min/med/max=%ld/%ld/%ld; vars:",
                 st, s->total, s->min_step_seconds, s->median_step_seconds, s->max_step_seconds);
             print_stage_variable_counts_inline(s);
             printf("\n");
        }

        if (active_stage >= 1 && active_stage <= 4 && sum.stage[active_stage].total > 0) {
            const CGEM_WQStageSummary *s = &sum.stage[active_stage];
            printf("  Active stage %d: last_target_step_seconds=%ld\n", active_stage, s->last_step_seconds);
        }

        if (WARMUP > 0 && active_stage >= 1 && active_stage <= 4) {
            printf("   Time base: target_step_seconds are interpreted as seconds since post-warmup start; sampling uses t_rel = (t - WARMUP) and does not record during warmup.\n");
        }
    }
}

static int g_detail_capture_log_remaining = 0;

static int get_param_index_stage1(const char *name)
{
    if (!name || !name[0]) return -1;
    if (strcasecmp(name, "Chezy1") == 0) return IDX1_CHEZY1;
    if (strcasecmp(name, "Chezy2") == 0) return IDX1_CHEZY2;
    if (strcasecmp(name, "LC1") == 0) return IDX1_LC1;
    if (strcasecmp(name, "LC2") == 0) return IDX1_LC2;
    if (strcasecmp(name, "Rs1") == 0) return IDX1_RS1;
    if (strcasecmp(name, "Rs2") == 0) return IDX1_RS2;
    if (strcasecmp(name, "C_VDB") == 0) return IDX1_C_VDB;
    if (strcasecmp(name, "D0_CORRECTION") == 0) return IDX1_D0_CORR;
    if (strcasecmp(name, "AMPL_CORRECTION") == 0) return IDX1_AMPL_CORR;
    return -1;
}

static double g1_lower_bounds[NUM_PARAMS_STAGE1];
static double g1_upper_bounds[NUM_PARAMS_STAGE1];
static double g1_initial_guess[NUM_PARAMS_STAGE1];
static double g1_best_params[NUM_PARAMS_STAGE1];
static double g1_best_score = DBL_MAX;
static ParamTransform g1_param_transform[NUM_PARAMS_STAGE1];

static double g2_lower_bounds[NUM_PARAMS_STAGE2];
static double g2_upper_bounds[NUM_PARAMS_STAGE2];
static double g2_initial_guess[NUM_PARAMS_STAGE2];
static double g2_best_params[NUM_PARAMS_STAGE2];
static double g2_best_score = DBL_MAX;
static ParamTransform g2_param_transform[NUM_PARAMS_STAGE2];

static double g3_lower_bounds[NUM_PARAMS_STAGE3];
static double g3_upper_bounds[NUM_PARAMS_STAGE3];
static double g3_initial_guess[NUM_PARAMS_STAGE3];
static double g3_best_params[NUM_PARAMS_STAGE3];
static double g3_best_score = DBL_MAX;
static ParamTransform g3_param_transform[NUM_PARAMS_STAGE3];

static double g4_lower_bounds[NUM_PARAMS_STAGE4];
static double g4_upper_bounds[NUM_PARAMS_STAGE4];
static double g4_initial_guess[NUM_PARAMS_STAGE4];
static double g4_best_params[NUM_PARAMS_STAGE4];
static double g4_best_score = DBL_MAX;
static ParamTransform g4_param_transform[NUM_PARAMS_STAGE4];
static int g_iteration_count = 0;
static int run_simulation_and_score(int stage, const double *params, double *score_out);
static double calibration_nlopt_objective_fn(unsigned n, const double *x, double *grad, void *data);

// Forward declaration (defined later)
static void ensure_output_directory(void);

// Forward declarations for stage application functions (defined later)
static void apply_parameters_stage1(const double *p);
static void apply_parameters_stage2(const double *p);
static void apply_parameters_stage3(const double *p);
static void apply_parameters_stage4(const double *p);

// Forward declarations for registry loaders and parameter name helpers (defined later)
static int load_parameter_registry_stage1(const char *path, double *lb, double *ub, double *x0, ParamTransform *xf);
static int load_parameter_registry_stage2(const char *path, double *lb, double *ub, double *x0, ParamTransform *xf);
static int load_parameter_registry_stage3(const char *path, double *lb, double *ub, double *x0, ParamTransform *xf);
static int load_parameter_registry_stage4(const char *path, double *lb, double *ub, double *x0, ParamTransform *xf);

static const char* get_param_name_stage1(int idx);
static const char* get_param_name_stage2(int idx);
static const char* get_param_name_stage3(int idx);
static const char* get_param_name_stage4(int idx);

typedef int (*LoadRegistryFn)(const char *path, double *lb, double *ub, double *x0, ParamTransform *xf);
typedef void (*ApplyParamsFn)(const double *p);
typedef const char* (*ParamNameFn)(int idx);

typedef struct {
    int stage;
    int num_params;
    double *lb;
    double *ub;
    double *x0;
    double *best_p;
    double *best_score;
    ParamTransform *xf;
    LoadRegistryFn load_registry;
    ApplyParamsFn apply_params;
    ParamNameFn param_name;
    const char *best_so_far_path;
} StageSpec;

static const StageSpec* get_stage_spec(int stage)
{
    switch (stage) {
        case 1: {
            static const StageSpec s = { 1, NUM_PARAMS_STAGE1,
                g1_lower_bounds, g1_upper_bounds, g1_initial_guess,
                g1_best_params, &g1_best_score, g1_param_transform,
                load_parameter_registry_stage1, apply_parameters_stage1,
                get_param_name_stage1,
                "OUT/Calibration/best_so_far_stage1.csv" };
            return &s;
        }
        case 2: {
            static const StageSpec s = { 2, NUM_PARAMS_STAGE2,
                g2_lower_bounds, g2_upper_bounds, g2_initial_guess,
                g2_best_params, &g2_best_score, g2_param_transform,
                load_parameter_registry_stage2, apply_parameters_stage2,
                get_param_name_stage2,
                "OUT/Calibration/best_so_far_stage2.csv" };
            return &s;
        }
        case 3: {
            static const StageSpec s = { 3, NUM_PARAMS_STAGE3,
                g3_lower_bounds, g3_upper_bounds, g3_initial_guess,
                g3_best_params, &g3_best_score, g3_param_transform,
                load_parameter_registry_stage3, apply_parameters_stage3,
                get_param_name_stage3,
                "OUT/Calibration/best_so_far_stage3.csv" };
            return &s;
        }
        case 4: {
            static const StageSpec s = { 4, NUM_PARAMS_STAGE4,
                g4_lower_bounds, g4_upper_bounds, g4_initial_guess,
                g4_best_params, &g4_best_score, g4_param_transform,
                load_parameter_registry_stage4, apply_parameters_stage4,
                get_param_name_stage4,
                "OUT/Calibration/best_so_far_stage4.csv" };
            return &s;
        }
        default:
            return NULL;
    }
}

static void save_best_so_far_row(const StageSpec *spec, double score, const double *params)
{
    if (!spec || !spec->best_so_far_path || !params) return;

    FILE *fp = fopen(spec->best_so_far_path, "a");
    if (!fp) return;

    fseek(fp, 0, SEEK_END);
    if (ftell(fp) == 0) {
        fprintf(fp, "iteration,score");
        for (int i = 0; i < spec->num_params; ++i) {
            const char *nm = spec->param_name ? spec->param_name(i) : NULL;
            if (nm && nm[0]) fprintf(fp, ",%s", nm);
            else fprintf(fp, ",param_%d", i);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "%d,%.10g", g_iteration_count, score);
    for (int i = 0; i < spec->num_params; ++i) {
        double val = params[i];
        if (spec->xf && spec->xf[i] == PARAM_TRANSFORM_LOG10) val = pow(10.0, val);
        fprintf(fp, ",%.10g", val);
    }
    fprintf(fp, "\n");

    fclose(fp);
}

static int get_param_index_stage2(const char *name)
{
    if (!name || !name[0]) return -1;
    if (strcasecmp(name, "Mero1") == 0) return IDX2_MERO1;
    if (strcasecmp(name, "tau_ero1") == 0) return IDX2_TAU_ERO1;
    if (strcasecmp(name, "tau_dep1") == 0) return IDX2_TAU_DEP1;
    if (strcasecmp(name, "Mero2") == 0) return IDX2_MERO2;
    if (strcasecmp(name, "tau_ero2") == 0) return IDX2_TAU_ERO2;
    if (strcasecmp(name, "tau_dep2") == 0) return IDX2_TAU_DEP2;
    if (strcasecmp(name, "Mero3") == 0) return IDX2_MERO3;
    if (strcasecmp(name, "tau_ero3") == 0) return IDX2_TAU_ERO3;
    if (strcasecmp(name, "tau_dep3") == 0) return IDX2_TAU_DEP3;
    if (strcasecmp(name, "Mero4") == 0) return IDX2_MERO4;
    if (strcasecmp(name, "tau_ero4") == 0) return IDX2_TAU_ERO4;
    if (strcasecmp(name, "tau_dep4") == 0) return IDX2_TAU_DEP4;
    if (strcasecmp(name, "P_ac") == 0) return IDX2_P_AC;
    if (strcasecmp(name, "K_ps") == 0) return IDX2_K_PS;
    if (strcasecmp(name, "k_ads") == 0) return IDX2_K_ADS;
    return -1;
}

static int get_param_index_stage3(const char *name)
{
    if (!name || !name[0]) return -1;
    if (strcasecmp(name, "scale_kox") == 0) return IDX3_SCALE_KOX;
    if (strcasecmp(name, "scale_knit") == 0) return IDX3_SCALE_KNIT;
    if (strcasecmp(name, "scale_kdenit") == 0) return IDX3_SCALE_KDENIT;
    if (strcasecmp(name, "scale_Pb_phy") == 0) return IDX3_SCALE_PB_PHY;
    if (strcasecmp(name, "scale_kmort_phy") == 0) return IDX3_SCALE_KMORT_PHY;
    if (strcasecmp(name, "scale_KD_Phy") == 0) return IDX3_SCALE_KD_PHY;
    if (strcasecmp(name, "scale_K_NH4_switch") == 0) return IDX3_SCALE_K_NH4_SWITCH;
    if (strcasecmp(name, "scale_kbg") == 0) return IDX3_SCALE_KBG;
    if (strcasecmp(name, "scale_kspm") == 0) return IDX3_SCALE_KSPM;
    if (strcasecmp(name, "scale_ws_phy") == 0) return IDX3_SCALE_WS_PHY;
    if (strcasecmp(name, "scale_alpha") == 0) return IDX3_SCALE_ALPHA;
    if (strcasecmp(name, "scale_KN") == 0) return IDX3_SCALE_KN;
    if (strcasecmp(name, "scale_KPO4") == 0) return IDX3_SCALE_KPO4;
    if (strcasecmp(name, "scale_KSi") == 0) return IDX3_SCALE_KSI;
    if (strcasecmp(name, "scale_KTOC") == 0) return IDX3_SCALE_KTOC;
    if (strcasecmp(name, "SOD_rate") == 0) return IDX3_SOD_RATE;
    if (strcasecmp(name, "kmort2_Phy1") == 0) return IDX3_KMORT2_PHY1;
    if (strcasecmp(name, "kcbod_fast") == 0) return IDX3_KCBOD_FAST;
    if (strcasecmp(name, "kCDOM") == 0) return IDX3_KCDOM;
    if (strcasecmp(name, "CN_toc") == 0) return IDX3_CN_TOC;
    return -1;
}


static int get_param_index_stage4(const char *name)
{
    if (!name || !name[0]) return -1;
    if (strcasecmp(name, "scale_piston_velocity") == 0) return IDX4_SCALE_PISTON_VEL;
    return -1;
}
// Current Value Getters (for default values if missing in CSV)

static double get_current_val_stage1(int idx) {
    switch (idx) {
        case IDX1_CHEZY1: return Chezy1;
        case IDX1_CHEZY2: return Chezy2;
        case IDX1_LC1: return LC1;
        case IDX1_LC2: return LC2;
        case IDX1_RS1: return Rs1;
        case IDX1_RS2: return Rs2;
        case IDX1_C_VDB: return C_VDB;
        case IDX1_D0_CORR: return D0_CORRECTION;
        case IDX1_AMPL_CORR: return AMPL_CORRECTION;
        default: return NAN;
    }
}

static double get_current_val_stage2(int idx) {
    switch (idx) {
        case IDX2_MERO1: return Mero1;
        case IDX2_TAU_ERO1: return tau_ero1;
        case IDX2_TAU_DEP1: return tau_dep1;
        case IDX2_MERO2: return Mero2;
        case IDX2_TAU_ERO2: return tau_ero2;
        case IDX2_TAU_DEP2: return tau_dep2;
        case IDX2_MERO3: return Mero3;
        case IDX2_TAU_ERO3: return tau_ero3;
        case IDX2_TAU_DEP3: return tau_dep3;
        case IDX2_MERO4: return Mero4;
        case IDX2_TAU_ERO4: return tau_ero4;
        case IDX2_TAU_DEP4: return tau_dep4;
        case IDX2_P_AC: return P_ac;
        case IDX2_K_PS: return K_ps;
        case IDX2_K_ADS: return k_ads;
        default: return NAN;
    }
}

static double get_current_val_stage3(int idx) {
    // Direct-value parameters: return current global value (set from params.txt)
    switch (idx) {
        case IDX3_SOD_RATE:    return SOD_rate;
        case IDX3_KMORT2_PHY1: return kmort2_Phy1;
        case IDX3_KCBOD_FAST:  return kcbod_fast;
        case IDX3_KCDOM:       return kCDOM;
        case IDX3_CN_TOC:      return CN_toc;
        default: break;
    }
    // Scale multipliers: default is 1.0 (applied on top of config values)
    return 1.0; 
}

static const char* get_param_name_stage1(int idx) {
    switch (idx) {
        case IDX1_CHEZY1: return "Chezy1";
        case IDX1_CHEZY2: return "Chezy2";
        case IDX1_LC1: return "LC1";
        case IDX1_LC2: return "LC2";
        case IDX1_RS1: return "Rs1";
        case IDX1_RS2: return "Rs2";
        case IDX1_C_VDB: return "C_VDB";
        case IDX1_D0_CORR: return "D0_CORRECTION";
        case IDX1_AMPL_CORR: return "AMPL_CORRECTION";
        default: return "UNKNOWN";
    }
}

static const char* get_param_name_stage2(int idx) {
    switch (idx) {
        case IDX2_MERO1: return "Mero1";
        case IDX2_TAU_ERO1: return "tau_ero1";
        case IDX2_TAU_DEP1: return "tau_dep1";
        case IDX2_MERO2: return "Mero2";
        case IDX2_TAU_ERO2: return "tau_ero2";
        case IDX2_TAU_DEP2: return "tau_dep2";
        case IDX2_MERO3: return "Mero3";
        case IDX2_TAU_ERO3: return "tau_ero3";
        case IDX2_TAU_DEP3: return "tau_dep3";
        case IDX2_MERO4: return "Mero4";
        case IDX2_TAU_ERO4: return "tau_ero4";
        case IDX2_TAU_DEP4: return "tau_dep4";
        case IDX2_P_AC: return "P_ac";
        case IDX2_K_PS: return "K_ps";
        case IDX2_K_ADS: return "k_ads";
        default: return "UNKNOWN";
    }
}

static const char* get_param_name_stage3(int idx) {
    switch (idx) {
        case IDX3_SCALE_KOX: return "scale_kox";
        case IDX3_SCALE_KNIT: return "scale_knit";
        case IDX3_SCALE_KDENIT: return "scale_kdenit";
        case IDX3_SCALE_PB_PHY: return "scale_Pb_phy";
        case IDX3_SCALE_KMORT_PHY: return "scale_kmort_phy";
        case IDX3_SCALE_KD_PHY: return "scale_KD_Phy";
        case IDX3_SCALE_K_NH4_SWITCH: return "scale_K_NH4_switch";
        case IDX3_SCALE_KBG: return "scale_kbg";
        case IDX3_SCALE_KSPM: return "scale_kspm";
        case IDX3_SCALE_WS_PHY: return "scale_ws_phy";
        case IDX3_SCALE_ALPHA: return "scale_alpha";
        case IDX3_SCALE_KN: return "scale_KN";
        case IDX3_SCALE_KPO4: return "scale_KPO4";
        case IDX3_SCALE_KSI: return "scale_KSi";
        case IDX3_SCALE_KTOC: return "scale_KTOC";
        case IDX3_SOD_RATE: return "SOD_rate";
        case IDX3_KMORT2_PHY1: return "kmort2_Phy1";
        case IDX3_KCBOD_FAST: return "kcbod_fast";
        case IDX3_KCDOM: return "kCDOM";
        case IDX3_CN_TOC: return "CN_toc";
        default: return "UNKNOWN";
    }
}


static const char* get_param_name_stage4(int idx) {
    switch (idx) {
        case IDX4_SCALE_PISTON_VEL: return "scale_piston_velocity";
        default: return "UNKNOWN";
    }
}

static double get_current_val_stage4(int idx) {
    (void)idx;
    return 1.0; 
}

// Generic Registry Loader

typedef int (*IndexLookupFn)(const char*);
typedef double (*CurrentValFn)(int);

static int load_registry_generic(const char *path, int stage, int num_params,
                               double *lb, double *ub, double *x0, ParamTransform *xf,
                               IndexLookupFn idx_lookup, CurrentValFn val_lookup)
{
    if (!path || !path[0]) return -1;
    if (!lb || !ub || !x0 || !xf || !idx_lookup) return -1;

    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "❌ ERROR: Cannot open %s\n", path);
        return -1;
    }

    bool *seen = (bool*)calloc((size_t)num_params, sizeof(bool));
    if (!seen) {
        fclose(fp);
        return -1;
    }

    char line[1024];
    int line_no = 0;
    while (fgets(line, sizeof(line), fp)) {
        line_no++;
        trim_trailing_whitespace(line);
        strip_utf8_bom(line);
        if (line[0] == '\0') continue;
        if (line[0] == '#') continue;

        // Skip header row
        if (strncmp(line, "stage,", 6) == 0) continue;

        char buf[1024];
        strncpy(buf, line, sizeof(buf) - 1);
        buf[sizeof(buf) - 1] = '\0';

        char *tok[10] = {0};
            if (!split_csv_row_n(buf, tok, 10)) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid row (expected exactly 10 columns).\n", path, line_no);
            free(seen);
            fclose(fp);
            return -1;
        }

        int row_stage = 0;
        if (!try_parse_int_strict(tok[0], &row_stage) || row_stage != stage) continue;

        const char *name = tok[1];
        if (!name || !name[0]) continue;

        int in_use = 1;
        if (!try_parse_enabled_flag(tok[8], &in_use)) {
            fprintf(stderr, "⚠️ WARNING: %s:%d invalid in_use '%s' for '%s'; defaulting to enabled\n",
                path, line_no, tok[8] ? tok[8] : "", name);
            in_use = 1;
        }

        int idx = idx_lookup(name);
        if (idx < 0 || idx >= num_params) {
            // Unknown parameter name for this stage: ignore (shared config file).
            continue;
        }

        if (seen[idx]) {
            fprintf(stderr, "❌ Duplicate parameter '%s' (stage %d)\n", name, stage);
            free(seen);
            fclose(fp);
            return -1;
        }

        double def_v = NAN, min_v = NAN, max_v = NAN;
        if (!try_parse_double_strict(tok[3], &def_v) || !isfinite(def_v)) def_v = NAN;
        if (!try_parse_double_strict(tok[4], &min_v) || !isfinite(min_v)) min_v = NAN;
        if (!try_parse_double_strict(tok[5], &max_v) || !isfinite(max_v)) max_v = NAN;

        if (!isfinite(min_v) || !isfinite(max_v) || !(max_v >= min_v)) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid bounds for '%s' (min/max).\n", path, line_no, name);
            free(seen);
            fclose(fp);
            return -1;
        }

        ParamTransform tf = PARAM_TRANSFORM_LINEAR;
        if (tok[7] && strcasecmp(tok[7], "log10") == 0) tf = PARAM_TRANSFORM_LOG10;

        // Determine initial guess: prefer CSV default_value (supports warm-start)
        double guess;
        if (isfinite(def_v)) {
            guess = def_v;
        } else if (val_lookup) {
            double cur = val_lookup(idx);
            guess = isfinite(cur) ? cur : min_v;
        } else {
            guess = min_v;
        }

        // Apply transform
        if (tf == PARAM_TRANSFORM_LOG10) {
            if (min_v <= 0 || max_v <= 0 || guess <= 0) {
                fprintf(stderr, "❌ ERROR: %s:%d log10 transform requires positive min/max/guess for '%s'\n", path, line_no, name);
                free(seen);
                fclose(fp);
                return -1;
            }
            min_v = log10(min_v);
            max_v = log10(max_v);
            guess = log10(guess);
        }

        // Clamp
        if (guess < min_v) guess = min_v;
        if (guess > max_v) guess = max_v;

        // If disabled, freeze at guess
        if (!in_use) {
            min_v = guess;
            max_v = guess;
        }

        lb[idx] = min_v;
        ub[idx] = max_v;
        x0[idx] = guess;
        xf[idx] = tf;
        seen[idx] = true;
    }

    fclose(fp);

    // Verify all required parameters found
    for (int i = 0; i < num_params; ++i) {
        if (!seen[i]) {
            fprintf(stderr, "❌ Error: Missing parameter index %d for Stage %d\n", i, stage);
            free(seen);
            return -1;
        }
    }

    free(seen);
    return 0;
}

static int load_parameter_registry_stage1(const char *path, double *lb, double *ub, double *x0, ParamTransform *xf) {
    return load_registry_generic(path, 1, NUM_PARAMS_STAGE1, lb, ub, x0, xf, get_param_index_stage1, get_current_val_stage1);
}

static int load_parameter_registry_stage2(const char *path, double *lb, double *ub, double *x0, ParamTransform *xf) {
    return load_registry_generic(path, 2, NUM_PARAMS_STAGE2, lb, ub, x0, xf, get_param_index_stage2, get_current_val_stage2);
}

static int load_parameter_registry_stage3(const char *path, double *lb, double *ub, double *x0, ParamTransform *xf) {
    return load_registry_generic(path, 3, NUM_PARAMS_STAGE3, lb, ub, x0, xf, get_param_index_stage3, get_current_val_stage3);
}

static int load_parameter_registry_stage4(const char *path, double *lb, double *ub, double *x0, ParamTransform *xf) {
    return load_registry_generic(path, 4, NUM_PARAMS_STAGE4, lb, ub, x0, xf, get_param_index_stage4, get_current_val_stage4);
}

static bool g_wq_needed[CHEM_COUNT];

// Scale factors (Stage 3 & 4)
static double scale_kox = 1.0;
static double scale_knit = 1.0;
static double scale_kdenit = 1.0;
static double scale_piston_velocity = 1.0;
static double scale_Pb_phy = 1.0;   // Unified max production multiplier
static double scale_kmort_phy = 1.0; // Unified mortality multiplier
static double scale_KD_Phy = 1.0;
static double scale_K_NH4_switch = 1.0;
static double scale_kbg = 1.0;   // Background light attenuation multiplier
static double scale_kspm = 1.0;  // SPM light attenuation multiplier
static double scale_ws_phy = 1.0;  // Phytoplankton settling velocity multiplier
static double scale_alpha = 1.0;  // Unified Alpha multiplier
static double scale_KN = 1.0;     // Nitrate half-saturation multiplier
static double scale_KPO4 = 1.0;   // Phosphate half-saturation multiplier
static double scale_KSi = 1.0;    // Silica half-saturation multiplier
static double scale_KTOC = 1.0;   // TOC half-saturation multiplier
static double calib_SOD_rate = NAN;            // Direct Stage-3 control
static double calib_kmort2_Phy1 = NAN;         // Direct Stage-3 control
static double calib_kcbod_fast = NAN;          // Direct Stage-3 control
static double calib_kCDOM = NAN;               // Direct Stage-3 CDOM absorption
static double calib_CN_toc = NAN;              // Direct Stage-3 C:N for TOC mineralization
// static double scale_C_Chla = 1.0; // C:Chla ratio multiplier

// Baseline (pre-calibration) biogeochemical constants to avoid compounding across
// iterations (especially when params.txt overrides are active).
static double g_base_Pb_phy1 = NAN;
static double g_base_Pb_phy2 = NAN;
static double g_base_kmort_phy1 = NAN;
static double g_base_kmort_phy2 = NAN;
static double g_base_kox = NAN;
static double g_base_knit = NAN;
static double g_base_kdenit = NAN;
static double g_base_KD_Phy = NAN;
static double g_base_K_NH4_switch = NAN;
static double g_base_piston_velocity_scale = NAN;
static double g_base_kbg = NAN;   // Baseline background light attenuation
static double g_base_kspm = NAN;  // Baseline SPM light attenuation
static double g_base_ws_phy = NAN; // Baseline phytoplankton settling velocity
static double g_base_alpha_phy1 = NAN; // Baseline alpha phy1
static double g_base_alpha_phy2 = NAN; // Baseline alpha phy2
static double g_base_KN_phy1 = NAN;    // Baseline KN
static double g_base_KN_phy2 = NAN;
static double g_base_KPO4_phy1 = NAN;  // Baseline KPO4
static double g_base_KPO4_phy2 = NAN;
static double g_base_KSi_phy1 = NAN;   // Baseline KSi
static double g_base_KSi_phy2 = NAN;
static double g_base_KTOC = NAN;       // Baseline KTOC
static double g_base_SOD_rate = NAN;   // Baseline SOD
static double g_base_kmort2_Phy1 = NAN;         // Baseline quadratic mortality closure
static double g_base_kcbod_fast = NAN;          // Baseline fast labile-C degradation
static double g_base_kCDOM = NAN;               // Baseline CDOM absorption
static double g_base_CN_toc = NAN;              // Baseline C:N ratio for TOC

static void capture_biogeo_baselines_after_params(void)
{
    // Ensure assignBiogeochemicalRateConstants() has been called at least once before
    // capturing baselines, otherwise many globals are still zero-initialized.
    g_base_Pb_phy1 = Pb[Phy1];
    g_base_Pb_phy2 = Pb[Phy2];
    g_base_kmort_phy1 = kmortality[Phy1];
    g_base_kmort_phy2 = kmortality[Phy2];
    g_base_kox = kox;
    g_base_knit = knit;
    g_base_kdenit = kdenit;
    g_base_KD_Phy = KD_Phy;
    g_base_K_NH4_switch = K_NH4_switch;
    g_base_piston_velocity_scale = piston_velocity_scale;
    g_base_kbg = kbg;
    g_base_kspm = kspm;
    g_base_ws_phy = ws_phy;
    g_base_alpha_phy1 = alpha[Phy1];
    g_base_alpha_phy2 = alpha[Phy2];
    g_base_KN_phy1 = KN[Phy1];
    g_base_KN_phy2 = KN[Phy2];
    g_base_KPO4_phy1 = KPO4[Phy1];
    g_base_KPO4_phy2 = KPO4[Phy2];
    g_base_KSi_phy1 = KSi[Phy1];
    g_base_KSi_phy2 = KSi[Phy2];
    g_base_KTOC = KTOC;
    g_base_SOD_rate = SOD_rate;
    g_base_kmort2_Phy1 = kmort2_Phy1;
    g_base_kcbod_fast = kcbod_fast;
    g_base_kCDOM = kCDOM;
    g_base_CN_toc = CN_toc;
}

static void restore_biogeo_baselines(void)
{
    if (isfinite(g_base_Pb_phy1)) Pb[Phy1] = g_base_Pb_phy1;
    if (isfinite(g_base_Pb_phy2)) Pb[Phy2] = g_base_Pb_phy2;
    if (isfinite(g_base_kmort_phy1)) kmortality[Phy1] = g_base_kmort_phy1;
    if (isfinite(g_base_kmort_phy2)) kmortality[Phy2] = g_base_kmort_phy2;
    if (isfinite(g_base_kox)) kox = g_base_kox;
    if (isfinite(g_base_knit)) knit = g_base_knit;
    if (isfinite(g_base_kdenit)) kdenit = g_base_kdenit;
    if (isfinite(g_base_KD_Phy)) KD_Phy = g_base_KD_Phy;
    if (isfinite(g_base_K_NH4_switch)) K_NH4_switch = g_base_K_NH4_switch;
    if (isfinite(g_base_piston_velocity_scale)) piston_velocity_scale = g_base_piston_velocity_scale;
    if (isfinite(g_base_kbg)) kbg = g_base_kbg;
    if (isfinite(g_base_kspm)) kspm = g_base_kspm;
    if (isfinite(g_base_ws_phy)) ws_phy = g_base_ws_phy;
    if (isfinite(g_base_alpha_phy1)) alpha[Phy1] = g_base_alpha_phy1;
    if (isfinite(g_base_alpha_phy2)) alpha[Phy2] = g_base_alpha_phy2;
    if (isfinite(g_base_KN_phy1)) KN[Phy1] = g_base_KN_phy1;
    if (isfinite(g_base_KN_phy2)) KN[Phy2] = g_base_KN_phy2;
    if (isfinite(g_base_KPO4_phy1)) KPO4[Phy1] = g_base_KPO4_phy1;
    if (isfinite(g_base_KPO4_phy2)) KPO4[Phy2] = g_base_KPO4_phy2;
    if (isfinite(g_base_KSi_phy1)) KSi[Phy1] = g_base_KSi_phy1;
    if (isfinite(g_base_KSi_phy2)) KSi[Phy2] = g_base_KSi_phy2;
    if (isfinite(g_base_KTOC)) KTOC = g_base_KTOC;
    if (isfinite(g_base_SOD_rate)) SOD_rate = g_base_SOD_rate;
    if (isfinite(g_base_kmort2_Phy1)) kmort2_Phy1 = g_base_kmort2_Phy1;
    if (isfinite(g_base_kcbod_fast)) kcbod_fast = g_base_kcbod_fast;
    if (isfinite(g_base_kCDOM)) kCDOM = g_base_kCDOM;
    if (isfinite(g_base_CN_toc)) CN_toc = g_base_CN_toc;
}

static void apply_biogeo_calibration_scales(int stage)
{
    // Start from baselines each iteration (prevents compounding when overrides are set).
    restore_biogeo_baselines();

    if (stage == 3) {
        if (scale_kox > 0.0 && isfinite(scale_kox)) kox *= scale_kox;
        if (scale_knit > 0.0 && isfinite(scale_knit)) knit *= scale_knit;
        if (scale_kdenit > 0.0 && isfinite(scale_kdenit)) kdenit *= scale_kdenit;
        if (scale_Pb_phy > 0.0 && isfinite(scale_Pb_phy)) {
            Pb[Phy1] *= scale_Pb_phy;
            Pb[Phy2] *= scale_Pb_phy;
        }
        if (scale_kmort_phy > 0.0 && isfinite(scale_kmort_phy)) {
            kmortality[Phy1] *= scale_kmort_phy;
            kmortality[Phy2] *= scale_kmort_phy;
        }

        if (scale_KD_Phy > 0.0 && isfinite(scale_KD_Phy)) KD_Phy *= scale_KD_Phy;
        if (scale_K_NH4_switch > 0.0 && isfinite(scale_K_NH4_switch)) K_NH4_switch *= scale_K_NH4_switch;

        // Light attenuation parameters (critical for phytoplankton dynamics)
        if (scale_kbg > 0.0 && isfinite(scale_kbg)) kbg *= scale_kbg;
        if (scale_kspm > 0.0 && isfinite(scale_kspm)) kspm *= scale_kspm;
        if (scale_ws_phy > 0.0 && isfinite(scale_ws_phy)) ws_phy *= scale_ws_phy;
        if (scale_alpha > 0.0 && isfinite(scale_alpha)) {
            alpha[Phy1] *= scale_alpha;
            alpha[Phy2] *= scale_alpha;
        }

        // Nutrient limitations
        if (scale_KN > 0.0 && isfinite(scale_KN)) {
            KN[Phy1] *= scale_KN;
            KN[Phy2] *= scale_KN;
        }
        if (scale_KPO4 > 0.0 && isfinite(scale_KPO4)) {
            KPO4[Phy1] *= scale_KPO4;
            KPO4[Phy2] *= scale_KPO4;
        }
        if (scale_KSi > 0.0 && isfinite(scale_KSi)) {
            KSi[Phy1] *= scale_KSi;
            KSi[Phy2] *= scale_KSi;
        }
        if (scale_KTOC > 0.0 && isfinite(scale_KTOC)) {
            KTOC *= scale_KTOC;
        }

        // Directly calibrated Stage-3 controls for O2 balance and spatial bloom confinement.
        if (isfinite(calib_SOD_rate) && calib_SOD_rate >= 0.0) {
            SOD_rate = calib_SOD_rate;
        }
        if (isfinite(calib_kmort2_Phy1) && calib_kmort2_Phy1 >= 0.0) {
            kmort2_Phy1 = calib_kmort2_Phy1;
        }
        if (isfinite(calib_kcbod_fast) && calib_kcbod_fast >= 0.0) {
            kcbod_fast = calib_kcbod_fast;
        }
        if (isfinite(calib_kCDOM) && calib_kCDOM >= 0.0) {
            kCDOM = calib_kCDOM;
        }
        if (isfinite(calib_CN_toc) && calib_CN_toc >= 1.0) {
            CN_toc = calib_CN_toc;
        }
    }


    if (stage == 4) {
        if (scale_piston_velocity > 0.0 && isfinite(scale_piston_velocity) && isfinite(g_base_piston_velocity_scale)) {
            piston_velocity_scale = g_base_piston_velocity_scale * scale_piston_velocity;
        }
    }
}
static void reset_wq_needed_flags(void)
{
    for (int s = 0; s < CHEM_COUNT; ++s) {
        g_wq_needed[s] = false;
    }
}

static void configure_wq_needed_from_targets(int current_stage)
{
    reset_wq_needed_flags();
    for (size_t i = 0; i < g_wq_seasonal.count; ++i) {
        const SeasonalWQTarget *p = &g_wq_seasonal.points[i];
        
        // Filter by stage: only enable variables required for the current calibration stage
        if (p->stage != current_stage) continue;

        if (p->kind == WQ_VAR_PHY_TOTAL) {
            g_wq_needed[Phy1] = true;
            g_wq_needed[Phy2] = true;
        } else {
            if (p->chem >= 0 && p->chem < CHEM_COUNT) {
                g_wq_needed[p->chem] = true;
            }
        }
    }
}

// Objective Functions for NLopt

static double calibration_nlopt_objective_fn(unsigned n, const double *x, double *grad, void *data)
{
    (void)n; (void)grad;
    const StageSpec *spec = (const StageSpec *)data;
    if (!spec || !x) return 1e9;

    double score = DBL_MAX;
    if (run_simulation_and_score(spec->stage, x, &score) != 0) return 1e9;

    if (isfinite(score) && spec->best_score && score < *(spec->best_score)) {
        *(spec->best_score) = score;
        if (spec->best_p) {
            for (int i = 0; i < spec->num_params; ++i) spec->best_p[i] = x[i];
        }

        if (debug_level >= DEBUG_LEVEL_ESSENTIAL) {
            printf("✅ [Stage%d] New best: %.6f\n", spec->stage, score);
        }

        save_best_so_far_row(spec, score, x);
    }

    return score;
}

// Fixed-schema loaders

static void free_tidal_targets(TidalTargetSet *s)
{
    if (!s) return;
    free(s->points);
    s->points = NULL;
    s->count = 0;
}

static void free_wq_seasonal_targets(SeasonalWQTargetSet *s)
{
    if (!s) return;
    free(s->points);
    s->points = NULL;
    s->count = 0;
}

static bool map_wq_variable_name(const char *name, WQVarKind *kind_out, Chem *chem_out, bool *is_chla_ugl_out)
{
    if (!name || !kind_out || !chem_out || !is_chla_ugl_out) return false;

    *is_chla_ugl_out = false;

    if (strcasecmp(name, "Chla_ugL") == 0 || strcasecmp(name, "Chla_ugl") == 0 || strcasecmp(name, "Chl-a_ugL") == 0) {
        *kind_out = WQ_VAR_PHY_TOTAL;
        *chem_out = Phy1;
        *is_chla_ugl_out = true;
        return true;
    }

    // Phytoplankton total (used as Chl-a proxy in many workflows).
    if (strcasecmp(name, "Phy") == 0 || strcasecmp(name, "Chla") == 0 || strcasecmp(name, "Chl-a") == 0 || strcasecmp(name, "ChlA") == 0) {
        *kind_out = WQ_VAR_PHY_TOTAL;
        *chem_out = Phy1;
        return true;
    }

    // Allow unified seasonal targets to include Sal/SPM as variables.
    if (strcasecmp(name, "Sal") == 0 || strcasecmp(name, "SAL") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = Sal; return true; }
    if (strcasecmp(name, "SPM") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = SPM; return true; }

    // WQ variables
    if (strcasecmp(name, "O2") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = O2; return true; }
    if (strcasecmp(name, "NH4") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = NH4; return true; }
    if (strcasecmp(name, "NO3") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = NO3; return true; }
    if (strcasecmp(name, "PO4") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = PO4; return true; }
    if (strcasecmp(name, "DSi") == 0 || strcasecmp(name, "Si") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = Si; return true; }
    if (strcasecmp(name, "TOC") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = TOC; return true; }
    if (strcasecmp(name, "DIC") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = DIC; return true; }
    if (strcasecmp(name, "AT") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = AT; return true; }
    if (strcasecmp(name, "pCO2") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = pCO2; return true; }
    if (strcasecmp(name, "pH") == 0 || strcasecmp(name, "PH") == 0) { *kind_out = WQ_VAR_CHEM; *chem_out = PH; return true; }

    return false;
}

static int load_wq_seasonal_targets_strict(const char *path, SeasonalWQTargetSet *out)
{
    if (!path || !path[0] || !out) return -1;
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "❌ ERROR: Cannot open WQ seasonal targets file: %s\n", path);
        return -1;
    }

    free_wq_seasonal_targets(out);

    char line[512];
    int line_no = 0;

    // Read header (strict)
    while (fgets(line, sizeof(line), fp)) {
        line_no++;
        trim_trailing_whitespace(line);
        strip_utf8_bom(line);
        if (line[0] == '\0') continue;
        if (line[0] == '#') continue;
        const char *expected9 = "target_step_seconds,cell_i0,cell_i1,alpha,stage,variable,target_value,scale,weight";
        if (strcmp(line, expected9) != 0) {
            fprintf(stderr, "❌ ERROR: Invalid seasonal targets header in %s\n   Expected: %s\n   Got:      %s\n",
                    path, expected9, line);
            fclose(fp);
            return -1;
        }
        break;
    }

    while (fgets(line, sizeof(line), fp)) {
        line_no++;
        trim_trailing_whitespace(line);
        strip_utf8_bom(line);
        if (line[0] == '\0') continue;
        if (line[0] == '#') continue;

        char buf[512];
        strncpy(buf, line, sizeof(buf) - 1);
        buf[sizeof(buf) - 1] = '\0';

        char *tok[9] = {0};
        if (!split_csv_row_n(buf, tok, 9)) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid row (expected exactly 9 columns).\n", path, line_no);
            fclose(fp);
            return -1;
        }

        SeasonalWQTarget p;
        memset(&p, 0, sizeof(p));
        p.recorded = false;
        p.model_value = NAN;
        p.stage = 0;
        p.is_chla_ugl = false;

        double step_d = NAN;
        if (!try_parse_double_strict(tok[0], &step_d) || step_d < 0.0) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid target_step_seconds: '%s'\n", path, line_no, tok[0]);
            fclose(fp);
            return -1;
        }
        if (step_d > (double)LONG_MAX) {
            fprintf(stderr, "❌ ERROR: %s:%d target_step_seconds too large\n", path, line_no);
            fclose(fp);
            return -1;
        }
        p.target_step_seconds = (long)llround(step_d);

        if (!try_parse_int_strict(tok[1], &p.cell_i0) || p.cell_i0 < 1) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid cell_i0: '%s'\n", path, line_no, tok[1]);
            fclose(fp);
            return -1;
        }
        if (!try_parse_int_strict(tok[2], &p.cell_i1) || p.cell_i1 < 1) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid cell_i1: '%s'\n", path, line_no, tok[2]);
            fclose(fp);
            return -1;
        }
        if (!try_parse_double_strict(tok[3], &p.alpha) || p.alpha < 0.0 || p.alpha > 1.0) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid alpha (expected 0..1): '%s'\n", path, line_no, tok[3]);
            fclose(fp);
            return -1;
        }

        const int col_stage = 4;
        const int col_var = 5;
        const int col_value = 6;
        const int col_scale = 7;
        const int col_weight = 8;

        if (!try_parse_int_strict(tok[col_stage], &p.stage) || p.stage < 1 || p.stage > 4) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid stage (expected 1..4): '%s'\n", path, line_no, tok[col_stage]);
            fclose(fp);
            return -1;
        }

        WQVarKind kind = WQ_VAR_CHEM;
        Chem chem = O2;
        bool is_chla_ugl = false;
        if (!map_wq_variable_name(tok[col_var], &kind, &chem, &is_chla_ugl)) {
            fprintf(stderr, "❌ ERROR: %s:%d unknown variable '%s'\n", path, line_no, tok[col_var]);
            fclose(fp);
            return -1;
        }
        p.kind = kind;
        p.chem = chem;
        p.is_chla_ugl = is_chla_ugl;

        if (!try_parse_double_strict(tok[col_value], &p.target_value) || p.target_value < 0.0) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid target_value (expected >=0): '%s'\n", path, line_no, tok[col_value]);
            fclose(fp);
            return -1;
        }
        if (!try_parse_double_strict(tok[col_scale], &p.scale) || !(p.scale > 0.0)) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid scale (expected >0): '%s'\n", path, line_no, tok[col_scale]);
            fclose(fp);
            return -1;
        }
        if (!try_parse_double_strict(tok[col_weight], &p.weight) || !(p.weight > 0.0)) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid weight (expected >0): '%s'\n", path, line_no, tok[col_weight]);
            fclose(fp);
            return -1;
        }

        if (p.cell_i1 < p.cell_i0) {
            fprintf(stderr, "❌ ERROR: %s:%d invalid cell indices (expected cell_i1>=cell_i0)\n", path, line_no);
            fclose(fp);
            return -1;
        }

        size_t new_count = out->count + 1;
        SeasonalWQTarget *tmp = (SeasonalWQTarget *)realloc(out->points, new_count * sizeof(SeasonalWQTarget));
        if (!tmp) {
            fprintf(stderr, "❌ ERROR: Out of memory while loading WQ seasonal targets\n");
            fclose(fp);
            return -1;
        }
        out->points = tmp;
        out->points[out->count] = p;
        out->count = new_count;
    }

    fclose(fp);

    if (out->count == 0) {
        fprintf(stderr, "❌ ERROR: WQ seasonal targets file contains no targets: %s\n", path);
        return -1;
    }
    return 0;
}

static int append_tidal_target_strict(TidalTargetSet *s, double location_km, double range_m)
{
    if (!s) return -1;

    // Strict: reject duplicate locations (no merging).
    const double tol = 1e-9;
    for (size_t i = 0; i < s->count; ++i) {
        if (fabs(s->points[i].location_km - location_km) <= tol) {
            fprintf(stderr, "❌ ERROR: Duplicate tidal target location_km=%.6f\n", location_km);
            return -1;
        }
    }

    size_t new_count = s->count + 1;
    TidalTarget *tmp = (TidalTarget *)realloc(s->points, new_count * sizeof(TidalTarget));
    if (!tmp) return -1;

    s->points = tmp;
    s->points[new_count - 1].location_km = location_km;
    s->points[new_count - 1].range_m = range_m;
    s->count = new_count;
    return 0;
}

static int load_tidal_targets_strict(const char *path, TidalTargetSet *out)
{
    if (!path || !out) return -1;

    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "❌ ERROR: Unable to open tidal targets file: %s\n", path);
        return -1;
    }

    char line[4096];
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        fprintf(stderr, "❌ ERROR: Tidal targets file is empty: %s\n", path);
        return -1;
    }
    trim_trailing_whitespace(line);
    strip_utf8_bom(line);

    // Strict: only accept canonical header.
    if (strcmp(line, "location_km,range_m") != 0) {
        fclose(fp);
        fprintf(stderr, "❌ ERROR: Invalid tidal targets header in %s\n", path);
        fprintf(stderr, "   Expected: location_km,range_m\n");
        fprintf(stderr, "   Got:      %s\n", line);
        return -1;
    }

    out->count = 0;
    out->points = NULL;

    while (fgets(line, sizeof(line), fp)) {
        trim_trailing_whitespace(line);
        if (line[0] == '\0' || line[0] == '#') continue;

        // Tokenize strict 2 columns and reject extras
        char *tokens[2] = {0};
        int token_count = 0;
        char *token = strtok(line, ",");
        while (token && token_count < 2) {
            tokens[token_count++] = token;
            token = strtok(NULL, ",");
        }
        if (token != NULL) {
            fclose(fp);
            fprintf(stderr, "❌ ERROR: Invalid tidal row (expected exactly 2 columns).\n");
            return -1;
        }
        if (token_count != 2) {
            fclose(fp);
            fprintf(stderr, "❌ ERROR: Invalid tidal row (expected exactly 2 columns).\n");
            return -1;
        }

        double location_km = NAN;
        double range_m = NAN;
        if (!try_parse_double_strict(tokens[0], &location_km) || location_km < 0.0) {
            fclose(fp);
            fprintf(stderr, "❌ ERROR: Invalid tidal location_km: '%s'\n", tokens[0]);
            return -1;
        }
        if (!try_parse_double_strict(tokens[1], &range_m) || range_m < 0.0) {
            fclose(fp);
            fprintf(stderr, "❌ ERROR: Invalid tidal range: '%s'\n", tokens[1]);
            return -1;
        }

        if (append_tidal_target_strict(out, location_km, range_m) != 0) {
            fclose(fp);
            fprintf(stderr, "❌ ERROR: Out of memory while loading tidal targets\n");
            return -1;
        }
    }

    fclose(fp);
    if (out->count == 0) {
        fprintf(stderr, "❌ ERROR: No tidal targets loaded from %s\n", path);
        return -1;
    }

    return 0;
}

// Interpolation / scoring

static double interpolate_tidal_range_profile_at_km(const double *range_profile, double location_km)
{
    if (!range_profile || M <= 1) {
        return (range_profile) ? range_profile[0] : NAN;
    }

    double dx_km = (double)DELXI / 1000.0;
    if (dx_km <= 0.0) {
        return range_profile[0];
    }

    double index = location_km / dx_km;
    if (index <= 0.0) {
        return range_profile[0];
    }
    double max_index = (double)(M - 1);
    if (index >= max_index) {
        return range_profile[M - 1];
    }

    int lower = (int)floor(index);
    int upper = lower + 1;
    double f = index - (double)lower;
    return range_profile[lower] + f * (range_profile[upper] - range_profile[lower]);
}

static double evaluate_tidal_rmse(void)
{
    if (g_tidal.count == 0 || M <= 0) {
        return DBL_MAX;
    }

    double *range_profile = (double *)malloc(sizeof(double) * (size_t)M);
    if (!range_profile) {
        return DBL_MAX;
    }

    for (int i = 1; i <= M; ++i) {
        double r = tidal_range_mean[i];
        if (!isfinite(r) || r < 0.0) r = 0.0;
        range_profile[i - 1] = r;
    }

    double sum_sq = 0.0;
    int count = 0;
    for (size_t i = 0; i < g_tidal.count; ++i) {
        double model = interpolate_tidal_range_profile_at_km(range_profile, g_tidal.points[i].location_km);
        double obs = g_tidal.points[i].range_m;
        double d = model - obs;
        sum_sq += d * d;
        count++;
    }

    free(range_profile);

    if (count <= 0) return DBL_MAX;
    return sqrt(sum_sq / (double)count);
}

static double safe_scaled_rmse(double rmse, double scale)
{
    if (!isfinite(rmse)) return 1e9;
    if (scale <= 0.0) return rmse;
    return rmse / scale;
}

typedef struct {
    int is_phy_total;
    Chem chem;
    CGEM_ScoreAgg agg;
} CGEM_ScoreBreakdownRow;

static int compare_score_breakdown_rows_desc_sum_w_sq(const void *a, const void *b)
{
    const CGEM_ScoreBreakdownRow *ra = (const CGEM_ScoreBreakdownRow *)a;
    const CGEM_ScoreBreakdownRow *rb = (const CGEM_ScoreBreakdownRow *)b;
    const double da = ra->agg.sum_w_sq;
    const double dbb = rb->agg.sum_w_sq;
    if (da < dbb) return 1;
    if (da > dbb) return -1;
    return 0;
}

// Simulation evaluation for a candidate parameter vector
static void apply_parameters_stage1(const double *p)
{
    Chezy1 = p[IDX1_CHEZY1];
    Chezy2 = p[IDX1_CHEZY2];
    LC1 = p[IDX1_LC1];
    LC2 = p[IDX1_LC2];
    Rs1 = p[IDX1_RS1];
    Rs2 = p[IDX1_RS2];
    C_VDB = p[IDX1_C_VDB];
    D0_CORRECTION = p[IDX1_D0_CORR];
    AMPL_CORRECTION = p[IDX1_AMPL_CORR];
}

static void apply_parameters_stage2(const double *p)
{
    Mero1 = p[IDX2_MERO1];
    tau_ero1 = p[IDX2_TAU_ERO1];
    tau_dep1 = p[IDX2_TAU_DEP1];
    Mero2 = p[IDX2_MERO2];
    tau_ero2 = p[IDX2_TAU_ERO2];
    tau_dep2 = p[IDX2_TAU_DEP2];
    Mero3 = p[IDX2_MERO3];
    tau_ero3 = p[IDX2_TAU_ERO3];
    tau_dep3 = p[IDX2_TAU_DEP3];
    Mero4 = p[IDX2_MERO4];
    tau_ero4 = p[IDX2_TAU_ERO4];
    tau_dep4 = p[IDX2_TAU_DEP4];
    P_ac = p[IDX2_P_AC];
    K_ps = p[IDX2_K_PS];
    k_ads = p[IDX2_K_ADS];
}

static void apply_parameters_stage3(const double *p)
{
    scale_kox = p[IDX3_SCALE_KOX];
    scale_knit = p[IDX3_SCALE_KNIT];
    scale_kdenit = p[IDX3_SCALE_KDENIT];
    scale_Pb_phy = p[IDX3_SCALE_PB_PHY];
    scale_kmort_phy = p[IDX3_SCALE_KMORT_PHY];
    scale_KD_Phy = p[IDX3_SCALE_KD_PHY];
    scale_K_NH4_switch = p[IDX3_SCALE_K_NH4_SWITCH];
    scale_kbg = p[IDX3_SCALE_KBG];     // Background light attenuation multiplier
    scale_kspm = p[IDX3_SCALE_KSPM];   // SPM light attenuation multiplier
    scale_ws_phy = p[IDX3_SCALE_WS_PHY]; // Phytoplankton settling velocity multiplier
    scale_alpha = p[IDX3_SCALE_ALPHA];
    scale_KN = p[IDX3_SCALE_KN];
    scale_KPO4 = p[IDX3_SCALE_KPO4];
    scale_KSi = p[IDX3_SCALE_KSI];
    scale_KTOC = p[IDX3_SCALE_KTOC];
    calib_SOD_rate = p[IDX3_SOD_RATE];
    calib_kmort2_Phy1 = p[IDX3_KMORT2_PHY1];
    calib_kcbod_fast = p[IDX3_KCBOD_FAST];
    calib_kCDOM = p[IDX3_KCDOM];
    calib_CN_toc = p[IDX3_CN_TOC];
    // NOTE: actual scaling is applied after assignBiogeochemicalRateConstants().
}


static void apply_parameters_stage4(const double *p)
{
    scale_piston_velocity = p[IDX4_SCALE_PISTON_VEL];
    // NOTE: actual scaling is applied after assignBiogeochemicalRateConstants().
}

static void ensure_out_calibration_dir(void)
{
    ensure_output_directory();
}

static void write_stage_overrides_file(int stage, const StageSpec *spec, const double *best_params_raw)
{
    if (!spec || !best_params_raw) return;

    ensure_out_calibration_dir();

    char path[256];
    snprintf(path, sizeof(path), "OUT/Calibration/params_overrides_stage%d.txt", stage);

    FILE *fp = fopen(path, "w");
    if (!fp) return;

    fprintf(fp, "# Auto-generated by CGEM calibration (stage %d)\n", stage);
    fprintf(fp, "# Format: key = value\n");

    // Decode transforms into physical parameter values.
    double decoded[32];
    for (int i = 0; i < spec->num_params; ++i) {
        double vraw = best_params_raw[i];
        if (spec->xf && spec->xf[i] == PARAM_TRANSFORM_LOG10) vraw = pow(10.0, vraw);
        decoded[i] = vraw;
    }

    if (stage == 1) {
        fprintf(fp, "Chezy1 = %.10g\n", decoded[IDX1_CHEZY1]);
        fprintf(fp, "Chezy2 = %.10g\n", decoded[IDX1_CHEZY2]);
        fprintf(fp, "LC1 = %.10g\n", decoded[IDX1_LC1]);
        fprintf(fp, "LC2 = %.10g\n", decoded[IDX1_LC2]);
        fprintf(fp, "Rs1 = %.10g\n", decoded[IDX1_RS1]);
        fprintf(fp, "Rs2 = %.10g\n", decoded[IDX1_RS2]);
        fprintf(fp, "C_VDB = %.10g\n", decoded[IDX1_C_VDB]);
        fprintf(fp, "D0_CORRECTION = %.10g\n", decoded[IDX1_D0_CORR]);
        fprintf(fp, "AMPL_CORRECTION = %.10g\n", decoded[IDX1_AMPL_CORR]);
    } else if (stage == 2) {
        fprintf(fp, "Mero1 = %.10g\n", decoded[IDX2_MERO1]);
        fprintf(fp, "tau_ero1 = %.10g\n", decoded[IDX2_TAU_ERO1]);
        fprintf(fp, "tau_dep1 = %.10g\n", decoded[IDX2_TAU_DEP1]);
        fprintf(fp, "Mero2 = %.10g\n", decoded[IDX2_MERO2]);
        fprintf(fp, "tau_ero2 = %.10g\n", decoded[IDX2_TAU_ERO2]);
        fprintf(fp, "tau_dep2 = %.10g\n", decoded[IDX2_TAU_DEP2]);
        fprintf(fp, "Mero3 = %.10g\n", decoded[IDX2_MERO3]);
        fprintf(fp, "tau_ero3 = %.10g\n", decoded[IDX2_TAU_ERO3]);
        fprintf(fp, "tau_dep3 = %.10g\n", decoded[IDX2_TAU_DEP3]);
        fprintf(fp, "Mero4 = %.10g\n", decoded[IDX2_MERO4]);
        fprintf(fp, "tau_ero4 = %.10g\n", decoded[IDX2_TAU_ERO4]);
        fprintf(fp, "tau_dep4 = %.10g\n", decoded[IDX2_TAU_DEP4]);
        fprintf(fp, "P_ac = %.10g\n", decoded[IDX2_P_AC]);
        fprintf(fp, "K_ps = %.10g\n", decoded[IDX2_K_PS]);
        fprintf(fp, "k_ads = %.10g\n", decoded[IDX2_K_ADS]);
    } else if (stage == 3) {
        // Convert multipliers into effective physical parameters using captured baselines.
        scale_kox = decoded[IDX3_SCALE_KOX];
        scale_knit = decoded[IDX3_SCALE_KNIT];
        scale_kdenit = decoded[IDX3_SCALE_KDENIT];
        scale_Pb_phy = decoded[IDX3_SCALE_PB_PHY];
        scale_kmort_phy = decoded[IDX3_SCALE_KMORT_PHY];
        scale_KD_Phy = decoded[IDX3_SCALE_KD_PHY];
        scale_K_NH4_switch = decoded[IDX3_SCALE_K_NH4_SWITCH];
        scale_kbg = decoded[IDX3_SCALE_KBG];
        scale_kspm = decoded[IDX3_SCALE_KSPM];
        scale_ws_phy = decoded[IDX3_SCALE_WS_PHY];
        scale_alpha = decoded[IDX3_SCALE_ALPHA];
        scale_KN = decoded[IDX3_SCALE_KN];
        scale_KPO4 = decoded[IDX3_SCALE_KPO4];
        scale_KSi = decoded[IDX3_SCALE_KSI];
        scale_KTOC = decoded[IDX3_SCALE_KTOC];
        calib_SOD_rate = decoded[IDX3_SOD_RATE];
        calib_kmort2_Phy1 = decoded[IDX3_KMORT2_PHY1];
        calib_kcbod_fast = decoded[IDX3_KCBOD_FAST];

        apply_biogeo_calibration_scales(3);

        fprintf(fp, "kox = %.10g\n", kox);
        fprintf(fp, "knit = %.10g\n", knit);
        fprintf(fp, "kdenit = %.10g\n", kdenit);
        fprintf(fp, "Pb_Phy1 = %.10g\n", Pb[Phy1]);
        fprintf(fp, "Pb_Phy2 = %.10g\n", Pb[Phy2]);
        fprintf(fp, "alpha_Phy1 = %.10g\n", alpha[Phy1]);
        fprintf(fp, "alpha_Phy2 = %.10g\n", alpha[Phy2]);
        fprintf(fp, "kmortality_Phy1 = %.10g\n", kmortality[Phy1]);
        fprintf(fp, "kmortality_Phy2 = %.10g\n", kmortality[Phy2]);
        fprintf(fp, "KD_Phy = %.10g\n", KD_Phy);
        fprintf(fp, "K_NH4_switch = %.10g\n", K_NH4_switch);
        fprintf(fp, "kbg = %.10g\n", kbg);
        fprintf(fp, "kspm = %.10g\n", kspm);
        fprintf(fp, "ws_phy = %.10g\n", ws_phy);
        fprintf(fp, "KN_Phy1 = %.10g\n", KN[Phy1]);
        fprintf(fp, "KN_Phy2 = %.10g\n", KN[Phy2]);
        fprintf(fp, "KPO4_Phy1 = %.10g\n", KPO4[Phy1]);
        fprintf(fp, "KPO4_Phy2 = %.10g\n", KPO4[Phy2]);
        fprintf(fp, "KSi_Phy1 = %.10g\n", KSi[Phy1]);
        fprintf(fp, "KSi_Phy2 = %.10g\n", KSi[Phy2]);
        fprintf(fp, "KTOC = %.10g\n", KTOC);
        fprintf(fp, "SOD_rate = %.10g\n", SOD_rate);
        fprintf(fp, "kmort2_Phy1 = %.10g\n", kmort2_Phy1);
        fprintf(fp, "kcbod_fast = %.10g\n", kcbod_fast);
        fprintf(fp, "kCDOM = %.10g\n", kCDOM);
        fprintf(fp, "CN_toc = %.10g\n", CN_toc);
    } else if (stage == 4) {
        scale_piston_velocity = decoded[IDX4_SCALE_PISTON_VEL];
        apply_biogeo_calibration_scales(4);
        fprintf(fp, "piston_velocity_scale = %.10g\n", piston_velocity_scale);
    }

    fclose(fp);
}

static void reset_accumulators(void)
{
    if (M <= 0) return;

    for (size_t i = 0; i < g_wq_seasonal.count; ++i) {
        g_wq_seasonal.points[i].recorded = false;
        g_wq_seasonal.points[i].model_value = NAN;
    }
}

static void capture_wq_seasonal_samples(long t)
{
    if (g_wq_seasonal.count == 0) return;

    // Do not sample during warmup. Seasonal targets are defined on a post-warmup
    // (analysis) clock, consistent with the forcing time base.
    if (t < WARMUP) return;

    const long t_rel = t - WARMUP;

    for (size_t i = 0; i < g_wq_seasonal.count; ++i) {
        SeasonalWQTarget *p = &g_wq_seasonal.points[i];
        if (p->recorded) continue;

        // target_step_seconds are seconds since post-warmup start.
        // Use the shared pure helper to lock the convention in tests.
        if (!cgem_should_capture_seasonal_target_at_sim_time(t, WARMUP, p->target_step_seconds)) continue;

        double v0 = NAN;
        double v1 = NAN;
        if (p->kind == WQ_VAR_PHY_TOTAL) {
            const double phy0 = v[Phy1].c[p->cell_i0] + v[Phy2].c[p->cell_i0];
            const double phy1 = v[Phy1].c[p->cell_i1] + v[Phy2].c[p->cell_i1];
            if (p->is_chla_ugl) {
                // Observation operator: Phy [mmolC/m^3] -> Chl-a [ug/L]
                // ug/L numerically equals mg/m^3.
                const double c_chl = (C_Chl_ratio > 0.0 && isfinite(C_Chl_ratio)) ? C_Chl_ratio : 15.0;
                const double to_chla = CHLA_MOLAR_MASS_C_MG_PER_MMOL / c_chl;
                v0 = phy0 * to_chla;
                v1 = phy1 * to_chla;
            } else {
                v0 = phy0;
                v1 = phy1;
            }
        } else {
            const int s = (int)p->chem;
            if (s < 0 || s >= CHEM_COUNT) {
                p->recorded = true;
                p->model_value = NAN;
                continue;
            }
            v0 = v[s].c[p->cell_i0];
            v1 = v[s].c[p->cell_i1];
        }

        if (isfinite(v0) && isfinite(v1)) {
            p->model_value = (1.0 - p->alpha) * v0 + p->alpha * v1;
        }
        p->recorded = true;

        if (debug_level >= DEBUG_LEVEL_DETAIL && g_detail_capture_log_remaining > 0 && p->stage == calibration_stage) {
            if (p->kind == WQ_VAR_PHY_TOTAL) {
                printf("[Calib Stage%d] capture Phy(total): t=%ld (t_rel=%ld) target_step_seconds=%ld cells=%d..%d alpha=%.3f model=%.6g\n",
                       calibration_stage, t, t_rel, p->target_step_seconds, p->cell_i0, p->cell_i1, p->alpha, p->model_value);
            } else {
                printf("[Calib Stage%d] capture %s: t=%ld (t_rel=%ld) target_step_seconds=%ld cells=%d..%d alpha=%.3f model=%.6g\n",
                       calibration_stage, chem_short_name(p->chem), t, t_rel, p->target_step_seconds, p->cell_i0, p->cell_i1, p->alpha, p->model_value);
            }
            g_detail_capture_log_remaining--;
        }
    }
}

static double evaluate_wq_seasonal_score(void)
{
    if (g_wq_seasonal.count == 0) return DBL_MAX;

    double sum_w = 0.0;
    double sum_w_sq = 0.0;

    const bool should_print_breakdown =
        (debug_level >= DEBUG_LEVEL_ESSENTIAL) &&
        (g_iteration_count == 1 || ((g_iteration_count % 50) == 0));

    CGEM_ScoreSample *samples = NULL;
    size_t sample_count = 0;
    if (should_print_breakdown) {
        // Worst-case: every target is eligible; allocate once.
        samples = (CGEM_ScoreSample *)malloc(g_wq_seasonal.count * sizeof(CGEM_ScoreSample));
        if (!samples) {
            // Printing is best-effort; never fail scoring due to OOM in diagnostics.
            sample_count = 0;
        }
    }

    for (size_t i = 0; i < g_wq_seasonal.count; ++i) {
        const SeasonalWQTarget *p = &g_wq_seasonal.points[i];
        if (!p->recorded || !isfinite(p->model_value)) continue;

        // Unified seasonal target file supports per-row stage.
        // Only score targets that match the current calibration stage.
        if (p->stage != calibration_stage) continue;

        const double d = (p->model_value - p->target_value) / p->scale;
        sum_w += p->weight;
        sum_w_sq += p->weight * d * d;

        if (samples) {
            CGEM_ScoreSample s;
            if (p->kind == WQ_VAR_PHY_TOTAL) {
                s.kind = CGEM_SCORE_SAMPLE_PHY_TOTAL;
                s.chem = Phy1;
            } else {
                s.kind = CGEM_SCORE_SAMPLE_CHEM;
                s.chem = p->chem;
            }
            s.d = d;
            s.weight = p->weight;
            samples[sample_count++] = s;
        }
    }

    if (!(sum_w > 0.0) || !isfinite(sum_w_sq)) {
        free(samples);
        return 1e9;
    }

    const double score = sqrt(sum_w_sq / sum_w);

    if (should_print_breakdown && samples && sample_count > 0) {
        CGEM_ScoreBreakdown bd;
        if (cgem_compute_score_breakdown(samples, sample_count, &bd) == 0) {
            CGEM_ScoreBreakdownRow rows[CHEM_COUNT + 1];
            int nrows = 0;

            for (int c = 0; c < CHEM_COUNT; ++c) {
                if (bd.by_chem[c].n > 0) {
                    rows[nrows].is_phy_total = 0;
                    rows[nrows].chem = (Chem)c;
                    rows[nrows].agg = bd.by_chem[c];
                    nrows++;
                }
            }
            if (bd.phy_total.n > 0) {
                rows[nrows].is_phy_total = 1;
                rows[nrows].chem = Phy1;
                rows[nrows].agg = bd.phy_total;
                nrows++;
            }

            // Sort by descending contribution to sum_w_sq.
            qsort(rows, (size_t)nrows, sizeof(CGEM_ScoreBreakdownRow), compare_score_breakdown_rows_desc_sum_w_sq);

            const int max_rows = (debug_level >= DEBUG_LEVEL_DETAIL) ? nrows : (nrows < 12 ? nrows : 12);

            printf("\n[Stage%d] Seasonal score breakdown (iter=%d, score=%.6f)\n", calibration_stage, g_iteration_count, score);
            printf("  %-12s %6s %12s %14s\n", "variable", "n", "rmse_scaled", "contrib_%");

            for (int irow = 0; irow < max_rows; ++irow) {
                const CGEM_ScoreBreakdownRow *r = &rows[irow];
                const char *name = r->is_phy_total ? "Phy(total)" : chem_short_name(r->chem);
                const double rmse = cgem_score_agg_rmse_scaled(&r->agg);
                const double pct = cgem_score_contribution_pct(r->agg.sum_w_sq, bd.total_sum_w_sq);
                printf("  %-12s %6zu %12.6f %13.2f\n", name, r->agg.n, rmse, pct);
            }

            if (max_rows < nrows && debug_level < DEBUG_LEVEL_DETAIL) {
                printf("  ... (%d more; set CGEM_DEBUG_LEVEL=%d for full list)\n", nrows - max_rows, DEBUG_LEVEL_DETAIL);
            }
        }
    }

    free(samples);
    return score;
}

static int run_simulation_and_score(int stage, const double *params, double *score_out)
{
    if (!params || !score_out) return -1;

    const StageSpec *spec = get_stage_spec(stage);
    if (!spec) return -1;

    g_iteration_count++;
    // Progress UI: use carriage return for interactive terminals, but also emit a
    // periodic newline so captured logs don't look like we "jump" from a small
    // iteration number straight to maxeval.
    printf("\r🔄 Iteration %-4d ", g_iteration_count);
    if (calibration_mode > 0 && (g_iteration_count % 50) == 0) {
        double best = (spec && spec->best_score) ? *(spec->best_score) : NAN;
        if (isfinite(best)) {
            printf("(best=%.6f)", best);
        }
        printf("\n");
    }
    fflush(stdout);

    reset_output_file_registry();
    // Calibration performance: do not delete OUT/* on every iteration.
    // We already suppress model outputs during calibration, and best-so-far rows are
    // appended under OUT/Calibration.
    if (calibration_mode == 0) {
        cleanupOldOutputFiles();
    }
    set_output_enabled(0);

    // Clear any fatal latch from a previous iteration.
    diagnostics_clear_fatal_request();

    // Backup biogeo state
    const int original_biogeo = enable_biogeochemical_reactions;
    
    // Configure Biology based on stage
    if (stage == 1) {
        enable_biogeochemical_reactions = BIOGEO_OFF; // Physics only
    } else if (stage == 2) {
        enable_biogeochemical_reactions = BIOGEO_OFF; // Stage 2 SPM-only: skip biogeo for speed
    } else {
        enable_biogeochemical_reactions = BIOGEO_FULL; // Needs full bio for Stage 3/4
    }

    // Reset multipliers defaults
    scale_kox = 1.0; scale_knit = 1.0; scale_kdenit = 1.0;
    scale_Pb_phy = 1.0; scale_kmort_phy = 1.0;
    scale_KD_Phy = 1.0; scale_K_NH4_switch = 1.0;
    scale_piston_velocity = 1.0;
    scale_kbg = 1.0; scale_kspm = 1.0; scale_ws_phy = 1.0; scale_alpha = 1.0;


    // Apply parameters BEFORE initializing physics
    const int n_params = spec->num_params;
    ParamTransform *transforms = spec->xf;

    // De-transform into physical parameter values
    double decoded[32];
    for (int i = 0; i < n_params; ++i) {
        decoded[i] = params[i];
        if (transforms && transforms[i] == PARAM_TRANSFORM_LOG10) decoded[i] = pow(10.0, params[i]);
    }

    if (spec->apply_params) spec->apply_params(decoded);

    diagnostics_init(debug_level);
    initializeHydrodynamics();
    initializeTransportVariables();
    initializeBiogeochemistry();

    if (diagnostics_fatal_requested()) {
        enable_biogeochemical_reactions = original_biogeo;
        set_output_enabled(1);
        if (score_out) *score_out = DBL_MAX;
        return -1;
    }
    assignBiogeochemicalRateConstants();
    // Apply Stage 3/4 scale factors on top of the freshly assigned constants.
    // This is the critical step that ensures calibration parameters actually affect
    // the dynamics being scored.
    apply_biogeo_calibration_scales(stage);

    // OPTIMIZATION: Disable transport for irrelevant variables based on stage
    // Save original env states
    int original_env[MAXV];
    for (int i=0; i<MAXV; ++i) original_env[i] = v[i].env;

    // Apply strict masking
    if (stage == 1) {
        // Stage 1: Tides + Salinity ONLY
        enable_suspended_sediment_dynamics = 0;
        enable_carbonate_diagnostics = 0; // explicit
        for (int i=0; i<MAXV; ++i) v[i].env = 0;
        v[Sal].env = 1;
        // Check if v[Sal] needs others? No.
    }
    else if (stage == 2) {
        // Stage 2: transport only Sal/SPM (PO4/PIP disabled for SPM-only speed)
        for (int i=0; i<MAXV; ++i) v[i].env = 0;
        v[Sal].env = 1;
        v[SPM].env = 1;
    }
    else if (stage == 3) {
        // Stage 3: eutrophication targets (typically nutrients, O2, phytoplankton).
        // Disable carbonate diagnostics unless explicitly targeted to avoid expensive
        // per-cell carbonate solves during optimization.
        enable_carbonate_diagnostics = 0;

        for (int i=0; i<MAXV; ++i) v[i].env = 0;

        // Always keep core transport/biogeo dependencies.
        v[Sal].env = 1; // mixing/transport baseline
        v[SPM].env = 1; // controls light attenuation
        v[Phy1].env = 1;
        v[Phy2].env = 1;
        v[NO3].env = 1;
        v[NH4].env = 1;
        v[PO4].env = 1;
        v[PIP].env = 1;
        v[Si].env = 1;
        v[O2].env = 1;
        v[TOC].env = 1;

        // Enable any additional variables requested by the seasonal target file.
        for (int s = 0; s < CHEM_COUNT; ++s) {
            if (g_wq_needed[s]) {
                v[s].env = 1;
            }
        }

        // If carbonate-related variables are targeted in Stage 3, turn diagnostics back on.
        if (g_wq_needed[pCO2] || g_wq_needed[PH] || g_wq_needed[CO2] || g_wq_needed[DIC] || g_wq_needed[AT]) {
            enable_carbonate_diagnostics = 1;
            v[pCO2].env = 1;
            v[PH].env = 1;
            v[CO2].env = 1;
            v[DIC].env = 1;
            v[AT].env = 1;
        }
    }
    else if (stage == 4) {
        // Stage 4: carbonate targets. Carbonate diagnostics and DIC/AT transport must be enabled.
        enable_carbonate_diagnostics = 1;

        for (int i=0; i<MAXV; ++i) v[i].env = 0;

        // Carbonate state + dependencies.
        v[Sal].env = 1;
        v[DIC].env = 1;
        v[AT].env = 1;
        v[pCO2].env = 1;
        v[PH].env = 1;
        v[CO2].env = 1;

        // Nutrients influence carbonate solver (NH4/NO3 terms).
        v[NH4].env = 1;
        v[NO3].env = 1;

        // Keep SPM/oxygen/TOC off by default unless targeted.
        for (int s = 0; s < CHEM_COUNT; ++s) {
            if (g_wq_needed[s]) {
                v[s].env = 1;
            }
        }
    }


    reset_accumulators();

    // Optional DETAIL diagnostics: show effective t_rel at the moment targets are captured.
    // To avoid log spam across long optimizations, only log a small number of captures
    // on the very first iteration.
    g_detail_capture_log_remaining = 0;
    if (debug_level >= DEBUG_LEVEL_DETAIL && g_iteration_count == 1) {
        g_detail_capture_log_remaining = 20;
    }

    long last_target_time = 0;
    for (size_t i = 0; i < g_wq_seasonal.count; ++i) {
        if (g_wq_seasonal.points[i].stage == stage) {
            if (g_wq_seasonal.points[i].target_step_seconds > last_target_time) {
                last_target_time = g_wq_seasonal.points[i].target_step_seconds;
            }
        }
    }

    int status = 0;
    for (long t = 0; t <= MAXT; t += DELTI) {
        // Optimization: Exit early once all seasonal targets for this stage are collected.
        if (t > WARMUP && (t - WARMUP) > last_target_time + DELTI) {
             break;
        }

        Hyd(t);
        bgboundary(t);
        Transport(t);
        if (diagnostics_fatal_requested()) {
            status = -1;
            break;
        }
        if (enable_suspended_sediment_dynamics) {
            updateSuspendedSediment(t);
        }

        if (diagnostics_fatal_requested()) {
            status = -1;
            break;
        }

        capture_wq_seasonal_samples(t);

        if (t > 0 && !check_numerical_stability(t)) {
            status = -1; break;
        }

        if (enable_biogeochemical_reactions) {
            Biogeo(t);
        }

        // In calibration mode, tracer invariant violations are latched (no exit).
        // Penalize this parameter set and move on.
        if (diagnostics_fatal_requested()) {
            status = -1;
            break;
        }
    }

    finalize_remaining_tidal_range();

    double score = DBL_MAX;

    if (status != 0) {
        score = DBL_MAX;
    } else {
        // Objective 2 only: unified long-form seasonal targets.
        if (stage == 1) {
            const double tidal = evaluate_tidal_rmse();
            const double sal_score = evaluate_wq_seasonal_score(); // stage=1, variable=Sal
            const double s_tidal = safe_scaled_rmse(tidal, 0.5);
            score = 0.4 * s_tidal + 0.6 * sal_score;
        } else {
            // Stage 2: seasonal SPM + PO4 (stage=2). Stage 3: seasonal WQ (stage=3). Stage 4: seasonal pCO2 (stage=4).
            score = evaluate_wq_seasonal_score();
        }
    }

    // Restore original env states
    for (int i=0; i<MAXV; ++i) v[i].env = original_env[i];
    enable_biogeochemical_reactions = original_biogeo;
    set_output_enabled(1);

    if (score_out) *score_out = score;
    return (status == 0 && isfinite(score)) ? 0 : -1;
}

// Public entry point
static void ensure_output_directory(void)
{
#ifdef _WIN32
    _mkdir("OUT");
    _mkdir("OUT/Calibration");
#else
    mkdir("OUT", 0777);
    mkdir("OUT/Calibration", 0777);
#endif
}

int run_hydrodynamic_calibration_from_params(void)
{
    fprintf(stderr, "[CALIB] run_hydrodynamic_calibration_from_params: begin\n");
    fflush(stderr);

    // Preflight: ensure the params path can be opened (helps localize hard crashes).
    {
        const char *pp = cgem_get_params_path();
        fprintf(stderr, "[CALIB] preflight params path: %s\n", pp ? pp : "(null)");
        fflush(stderr);
        FILE *fp = fopen(pp, "r");
        if (!fp) {
            fprintf(stderr, "[CALIB] preflight fopen failed for params path: %s\n", pp ? pp : "(null)");
            fflush(stderr);
        } else {
            fclose(fp);
            fprintf(stderr, "[CALIB] preflight fopen OK: %s\n", pp ? pp : "(null)");
            fflush(stderr);
        }
    }

    // Load full configuration
    read_parameters(cgem_get_params_path());

    fprintf(stderr, "[CALIB] read_parameters: OK (MAXT=%ld WARMUP=%ld DELTI=%d)\n", MAXT, WARMUP, DELTI);
    fflush(stderr);

    // Establish baseline biogeochemical constants from params.txt + defaults.
    // This ensures Stage 3/4 multipliers are applied to the correct base values and
    // do not compound across iterations.
    assignBiogeochemicalRateConstants();
    capture_biogeo_baselines_after_params();

    fprintf(stderr, "[CALIB] assignBiogeochemicalRateConstants + capture baselines: OK\n");
    fflush(stderr);

    if (calibration_stage < 1 || calibration_stage > 4) {
        fprintf(stderr, "❌ ERROR: calibration_stage=%d not supported (1-4 only).\n", calibration_stage);
        return -1;
    }

    if (calibration_mode == 0) {
        printf("⚠️ calibration_mode=0 disables optimization; skipping calibration.\n");
        return 0;
    }

    ensure_output_directory();

    fprintf(stderr, "[CALIB] ensured OUT/Calibration\n");
    fflush(stderr);

    // Load Targets
    free_tidal_targets(&g_tidal);
    free_wq_seasonal_targets(&g_wq_seasonal);

    // Unified seasonal objective: one long-form seasonal target file with per-row stage.
    // Stage 1 still uses tidal range targets.
    if (calibration_stage == 1) {
        if (load_tidal_targets_strict(calibration_tidal_targets_file, &g_tidal) != 0) return -1;
    }

    fprintf(stderr, "[CALIB] loaded tidal targets (stage1=%s)\n", (calibration_stage == 1) ? "yes" : "no");
    fflush(stderr);

    if (load_wq_seasonal_targets_strict(calibration_wq_seasonal_targets_file, &g_wq_seasonal) != 0) return -1;
    configure_wq_needed_from_targets(calibration_stage);

    fprintf(stderr, "[CALIB] loaded seasonal targets: n=%zu\n", g_wq_seasonal.count);
    fflush(stderr);

    // Phase 1 diagnostics: summarize seasonal targets and flag possible time-base mismatch.
    log_loaded_seasonal_targets_summary(calibration_stage);

    fprintf(stderr, "[CALIB] logged seasonal target summary\n");
    fflush(stderr);

    // Load Model Inputs
    readRiverbedProfile(riverbed_depth, riverbed_profile_file);
    allocateCumulativeDischargeMemory();
    readGlobalConfigSettings("INPUT/config_input.txt");
    readBoundaryData("INPUT/config_input.txt");
    if (tributaryEnabled) readTributaryData("INPUT/config_input.txt");
    bgboundary(0);

    fprintf(stderr, "[CALIB] loaded riverbed/config/boundaries (tributaryEnabled=%d)\n", tributaryEnabled);
    fflush(stderr);

    // Initialize Parameters based on Stage
    const StageSpec *spec = get_stage_spec(calibration_stage);
    if (!spec) return -1;

    const int num_params = spec->num_params;
    double *lb = spec->lb;
    double *ub = spec->ub;
    double *x0 = spec->x0;
    double *best_p = spec->best_p;
    double *best_score_ptr = spec->best_score;
    ParamTransform *xf = spec->xf;

    if (spec->load_registry && spec->load_registry(calibration_parameters_file, lb, ub, x0, xf) != 0) return -1;

    fprintf(stderr, "[CALIB] loaded calibration parameter registry for stage %d (num_params=%d)\n", calibration_stage, num_params);
    {
        int free_dims = 0;
        int fixed_dims = 0;
        for (int i = 0; i < num_params; ++i) {
            if (ub[i] > lb[i]) free_dims++;
            else fixed_dims++;
        }
        fprintf(stderr, "[CALIB] parameter DOF summary: free=%d fixed=%d\n", free_dims, fixed_dims);
    }
    fflush(stderr);

    // Init tracking
    if (best_score_ptr) *best_score_ptr = DBL_MAX;
    for (int i = 0; i < num_params; ++i) best_p[i] = x0[i];

    // Optimizer Setup
    nlopt_algorithm alg = NLOPT_LN_BOBYQA;
    const char* alg_name = "BOBYQA (Local)";

    switch(calibration_mode) {
        case 1: alg = NLOPT_LN_BOBYQA; alg_name = "BOBYQA (Local)"; break;
        case 2: alg = NLOPT_LN_NELDERMEAD; alg_name = "Nelder-Mead (Local)"; break;
        case 3: alg = NLOPT_LN_SBPLX; alg_name = "SBPLX (Local)"; break;
        case 4: alg = NLOPT_LN_COBYLA; alg_name = "COBYLA (Local)"; break;
        case 5: alg = NLOPT_GN_DIRECT_L; alg_name = "DIRECT-L (Global)"; break;
        case 6: alg = NLOPT_GN_CRS2_LM; alg_name = "CRS2-LM (Global)"; break;
        case 7: alg = NLOPT_GN_ISRES; alg_name = "ISRES (Global)"; break;
        default: alg = NLOPT_LN_BOBYQA; alg_name = "BOBYQA (Default/Local)"; break;
    }

    printf("\n🚀 Starting Calibration Stage %d using %s\n", calibration_stage, alg_name);
    g_iteration_count = 0;

    fprintf(stderr, "[CALIB] creating NLopt optimizer...\n");
    fflush(stderr);

    nlopt_opt opt = nlopt_create(alg, (unsigned)num_params);
    if (!opt) return -1;

    fprintf(stderr, "[CALIB] NLopt optimizer created\n");
    fflush(stderr);

    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, calibration_nlopt_objective_fn, (void*)spec);
    nlopt_set_maxeval(opt, (calibration_max_iterations > 0) ? calibration_max_iterations : 200);
    // For global modes, avoid tolerance-based early stopping so runs honor maxeval.
    if (calibration_mode == 5 || calibration_mode == 6 || calibration_mode == 7) {
        nlopt_set_xtol_rel(opt, 0.0);
        nlopt_set_ftol_rel(opt, 0.0);
    } else {
        nlopt_set_xtol_rel(opt, 1e-8);
        nlopt_set_ftol_rel(opt, 1e-8);
    }

    // Initial X
    double *x = (double*)malloc(sizeof(double)*num_params);
    if (!x) {
        nlopt_destroy(opt);
        return -1;
    }
    for(int i=0; i<num_params; ++i) x[i] = x0[i];

    // Perturb the initial guess slightly (5% of range) to avoid degenerate SBPLX/BOBYQA
    // simplex at the default all-ones position. Uses a deterministic hash-based jitter.
    for (int i = 0; i < num_params; ++i) {
        double range = ub[i] - lb[i];
        if (range > 0.0) {
            // Deterministic small perturbation: i-dependent direction
            double jitter = range * 0.05 * ((i % 2 == 0) ? 1.0 : -1.0);
            double perturbed = x[i] + jitter;
            if (perturbed < lb[i]) perturbed = lb[i] + range * 0.01;
            if (perturbed > ub[i]) perturbed = ub[i] - range * 0.01;
            x[i] = perturbed;
        }
    }

    double minf = DBL_MAX;
    fprintf(stderr, "[CALIB] calling nlopt_optimize (maxeval=%d, num_params=%d)...\n",
            (calibration_max_iterations > 0) ? calibration_max_iterations : 200, num_params);
    fflush(stderr);
    nlopt_result r = nlopt_optimize(opt, x, &minf);
    fprintf(stderr, "[CALIB] nlopt_optimize returned r=%d minf=%g evals=%d\n", (int)r, minf, g_iteration_count);
    fflush(stderr);

    nlopt_destroy(opt);
    free(x);

    // Hybrid improvement: global optimizers explore well but often plateau.
    // Follow up with a local refinement from the best point found so far.
    const bool is_global = (calibration_mode == 5 || calibration_mode == 6 || calibration_mode == 7);
    if (is_global && best_score_ptr && isfinite(*best_score_ptr) && *best_score_ptr < DBL_MAX) {
        const int refine_maxeval = (calibration_max_iterations > 0) ? (int)fmax(50.0, fmin(300.0, (double)calibration_max_iterations * 0.2)) : 200;
        printf("\n🛠️  Refining best solution with BOBYQA (maxeval=%d)...\n", refine_maxeval);

        nlopt_opt opt2 = nlopt_create(NLOPT_LN_BOBYQA, (unsigned)num_params);
        if (opt2) {
            nlopt_set_lower_bounds(opt2, lb);
            nlopt_set_upper_bounds(opt2, ub);
            nlopt_set_min_objective(opt2, calibration_nlopt_objective_fn, (void*)spec);
            nlopt_set_maxeval(opt2, refine_maxeval);
            nlopt_set_xtol_rel(opt2, 0);

            double *x2 = (double*)malloc(sizeof(double)*num_params);
            if (x2) {
                for (int i = 0; i < num_params; ++i) x2[i] = best_p[i];
                double minf2 = DBL_MAX;
                nlopt_result r2 = nlopt_optimize(opt2, x2, &minf2);
                if (r2 < 0) {
                    fprintf(stderr, "⚠️ Local refinement terminated with status %d\n", r2);
                }
                free(x2);
            }

            nlopt_destroy(opt2);
        }
    }

    if (r < 0) {
        // NLopt can return a negative status (e.g., -4 roundoff-limited) even when a
        // good solution was found. Do not drop the best parameters in that case.
        fprintf(stderr, "⚠️ NLopt terminated with status %d\n", r);
    }

    if (*best_score_ptr < DBL_MAX) {
        printf("\n✅ [Stage %d] Optimization Complete. Best Score: %.6f\n", calibration_stage, *best_score_ptr);
        // Dump best params (with names) and emit a params.txt override file.
        FILE *fp = fopen("OUT/Calibration/calibrated_params.csv", "w");
        if(fp) {
            fprintf(fp, "stage,index,param_name,value\n");
            for(int i=0; i<num_params; ++i) {
                double val = best_p[i];
                if (xf[i] == PARAM_TRANSFORM_LOG10) val = pow(10.0, best_p[i]);
                const char *nm = (spec && spec->param_name) ? spec->param_name(i) : NULL;
                if (!nm) nm = "param";
                fprintf(fp, "%d,%d,%s,%.10g\n", calibration_stage, i, nm, val);
            }
            fclose(fp);
        }

        write_stage_overrides_file(calibration_stage, spec, best_p);
    }

    return (*best_score_ptr < DBL_MAX) ? 0 : -1;
}
