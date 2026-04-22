// Windows/MSVC NLopt runtime shim
//
// Purpose
// -------
// The repository vendors an NLopt runtime DLL under external/nlopt/build.
// That DLL may come with a MinGW import library (libnlopt.dll.a), which MSVC
// cannot link against. To keep calibration runnable on Windows with MSVC,
// we provide small wrapper definitions for the subset of NLopt API used by
// src/calibration.c, and resolve them dynamically from libnlopt.dll at runtime.
//
// This file is only compiled when CGEM_ENABLE_CALIBRATION=ON and MSVC is used.

#if defined(_WIN32)

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>

#include <nlopt.h>

// -----------------------------------------------------------------------------
// Dynamic loader
// -----------------------------------------------------------------------------

static HMODULE g_nlopt_dll = NULL;

static void *load_symbol(const char *name)
{
    if (!g_nlopt_dll || !name) {
        return NULL;
    }

    FARPROC p = GetProcAddress(g_nlopt_dll, name);
    return (void *)p;
}

static int ensure_nlopt_loaded(void)
{
    if (g_nlopt_dll) {
        return 1;
    }

    // Optional explicit path override
    const char *env_path = getenv("CGEM_NLOPT_DLL");
    if (env_path && env_path[0] != '\0') {
        g_nlopt_dll = LoadLibraryA(env_path);
        if (g_nlopt_dll) {
            return 1;
        }
    }

    // Try standard search order first (DLL next to exe, PATH, system dirs)
    g_nlopt_dll = LoadLibraryA("libnlopt.dll");
    if (g_nlopt_dll) {
        return 1;
    }

    // Try repo-relative default location when running from repository root
    g_nlopt_dll = LoadLibraryA("external\\nlopt\\build\\libnlopt.dll");
    if (g_nlopt_dll) {
        return 1;
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "================================================================================\n");
    fprintf(stderr, "❌ NLopt runtime not found (calibration cannot run)\n");
    fprintf(stderr, "================================================================================\n");
    fprintf(stderr, "Tried to load: libnlopt.dll (PATH / next-to-exe)\n");
    fprintf(stderr, "Also tried: external/nlopt/build/libnlopt.dll (repo-relative)\n");
    fprintf(stderr, "Option: set environment variable CGEM_NLOPT_DLL to the full path of libnlopt.dll\n");
    fprintf(stderr, "================================================================================\n");
    fprintf(stderr, "\n");

    return 0;
}

// -----------------------------------------------------------------------------
// NLopt API subset (wrappers)
// -----------------------------------------------------------------------------

typedef nlopt_opt (*fn_nlopt_create)(nlopt_algorithm algorithm, unsigned n);
typedef void (*fn_nlopt_destroy)(nlopt_opt opt);
typedef nlopt_result (*fn_nlopt_set_lower_bounds)(nlopt_opt opt, const double *lb);
typedef nlopt_result (*fn_nlopt_set_upper_bounds)(nlopt_opt opt, const double *ub);
typedef nlopt_result (*fn_nlopt_set_min_objective)(nlopt_opt opt, nlopt_func f, void *f_data);
typedef nlopt_result (*fn_nlopt_set_maxeval)(nlopt_opt opt, int maxeval);
typedef nlopt_result (*fn_nlopt_set_population)(nlopt_opt opt, unsigned pop);
typedef nlopt_result (*fn_nlopt_set_xtol_rel)(nlopt_opt opt, double tol);
typedef nlopt_result (*fn_nlopt_set_ftol_rel)(nlopt_opt opt, double tol);
typedef nlopt_result (*fn_nlopt_optimize)(nlopt_opt opt, double *x, double *opt_f);

static fn_nlopt_create p_nlopt_create = NULL;
static fn_nlopt_destroy p_nlopt_destroy = NULL;
static fn_nlopt_set_lower_bounds p_nlopt_set_lower_bounds = NULL;
static fn_nlopt_set_upper_bounds p_nlopt_set_upper_bounds = NULL;
static fn_nlopt_set_min_objective p_nlopt_set_min_objective = NULL;
static fn_nlopt_set_maxeval p_nlopt_set_maxeval = NULL;
static fn_nlopt_set_population p_nlopt_set_population = NULL;
static fn_nlopt_set_xtol_rel p_nlopt_set_xtol_rel = NULL;
static fn_nlopt_set_ftol_rel p_nlopt_set_ftol_rel = NULL;
static fn_nlopt_optimize p_nlopt_optimize = NULL;

static int ensure_symbols_loaded(void)
{
    if (p_nlopt_create) {
        return 1;
    }

    if (!ensure_nlopt_loaded()) {
        return 0;
    }

    p_nlopt_create = (fn_nlopt_create)load_symbol("nlopt_create");
    p_nlopt_destroy = (fn_nlopt_destroy)load_symbol("nlopt_destroy");
    p_nlopt_set_lower_bounds = (fn_nlopt_set_lower_bounds)load_symbol("nlopt_set_lower_bounds");
    p_nlopt_set_upper_bounds = (fn_nlopt_set_upper_bounds)load_symbol("nlopt_set_upper_bounds");
    p_nlopt_set_min_objective = (fn_nlopt_set_min_objective)load_symbol("nlopt_set_min_objective");
    p_nlopt_set_maxeval = (fn_nlopt_set_maxeval)load_symbol("nlopt_set_maxeval");
    p_nlopt_set_population = (fn_nlopt_set_population)load_symbol("nlopt_set_population");
    p_nlopt_set_xtol_rel = (fn_nlopt_set_xtol_rel)load_symbol("nlopt_set_xtol_rel");
    p_nlopt_set_ftol_rel = (fn_nlopt_set_ftol_rel)load_symbol("nlopt_set_ftol_rel");
    p_nlopt_optimize = (fn_nlopt_optimize)load_symbol("nlopt_optimize");

    if (!p_nlopt_create || !p_nlopt_destroy || !p_nlopt_set_lower_bounds || !p_nlopt_set_upper_bounds ||
        !p_nlopt_set_min_objective || !p_nlopt_set_maxeval || !p_nlopt_set_population || !p_nlopt_set_xtol_rel ||
        !p_nlopt_set_ftol_rel || !p_nlopt_optimize) {
        fprintf(stderr, "❌ NLopt DLL loaded, but required symbols are missing. Calibration cannot run.\n");
        return 0;
    }

    return 1;
}

nlopt_opt nlopt_create(nlopt_algorithm algorithm, unsigned n)
{
    if (!ensure_symbols_loaded()) {
        return NULL;
    }
    return p_nlopt_create(algorithm, n);
}

void nlopt_destroy(nlopt_opt opt)
{
    if (!ensure_symbols_loaded()) {
        return;
    }
    p_nlopt_destroy(opt);
}

nlopt_result nlopt_set_lower_bounds(nlopt_opt opt, const double *lb)
{
    if (!ensure_symbols_loaded()) {
        return NLOPT_FAILURE;
    }
    return p_nlopt_set_lower_bounds(opt, lb);
}

nlopt_result nlopt_set_upper_bounds(nlopt_opt opt, const double *ub)
{
    if (!ensure_symbols_loaded()) {
        return NLOPT_FAILURE;
    }
    return p_nlopt_set_upper_bounds(opt, ub);
}

nlopt_result nlopt_set_min_objective(nlopt_opt opt, nlopt_func f, void *f_data)
{
    if (!ensure_symbols_loaded()) {
        return NLOPT_FAILURE;
    }
    return p_nlopt_set_min_objective(opt, f, f_data);
}

nlopt_result nlopt_set_maxeval(nlopt_opt opt, int maxeval)
{
    if (!ensure_symbols_loaded()) {
        return NLOPT_FAILURE;
    }
    return p_nlopt_set_maxeval(opt, maxeval);
}

nlopt_result nlopt_set_population(nlopt_opt opt, unsigned pop)
{
    if (!ensure_symbols_loaded()) {
        return NLOPT_FAILURE;
    }
    return p_nlopt_set_population(opt, pop);
}

nlopt_result nlopt_set_xtol_rel(nlopt_opt opt, double tol)
{
    if (!ensure_symbols_loaded()) {
        return NLOPT_FAILURE;
    }
    return p_nlopt_set_xtol_rel(opt, tol);
}

nlopt_result nlopt_set_ftol_rel(nlopt_opt opt, double tol)
{
    if (!ensure_symbols_loaded()) {
        return NLOPT_FAILURE;
    }
    return p_nlopt_set_ftol_rel(opt, tol);
}

nlopt_result nlopt_optimize(nlopt_opt opt, double *x, double *opt_f)
{
    if (!ensure_symbols_loaded()) {
        return NLOPT_FAILURE;
    }
    return p_nlopt_optimize(opt, x, opt_f);
}

#endif // _WIN32
