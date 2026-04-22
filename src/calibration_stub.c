/**
 * @file calibration_stub.c
 * @brief Stub implementation of calibration module without NLopt dependency
 * @author Nguyen Truong An
 * @date 2026-02-01
 *
 * This stub allows the core C-GEM model to build and run without requiring
 * NLopt or calibration support. Used when CGEM_ENABLE_CALIBRATION=0.
 *
 * SCIENTIFIC GUARDRAILS: This is a build/portability change only.
 * Does not affect hydrodynamic/transport physics implementation.
 */

#include <stdio.h>
#include <stdlib.h>
#include "calibration.h"
#include "define.h"

/**
 * @brief Stub for hydrodynamic calibration - prints error and returns failure code
 * @return -1 (calibration not available)
 *
 * This stub is called when calibration_mode != 0 but calibration support
 * was not compiled in (CGEM_ENABLE_CALIBRATION=0).
 */
int run_hydrodynamic_calibration_from_params(void) {
    fprintf(stderr, "\n");
    fprintf(stderr, "================================================================================\n");
    fprintf(stderr, "❌ CALIBRATION NOT AVAILABLE\n");
    fprintf(stderr, "================================================================================\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "This build was compiled WITHOUT calibration support.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "To enable calibration:\n");
    fprintf(stderr, "  1. Install NLopt library (https://nlopt.readthedocs.io/)\n");
    fprintf(stderr, "  2. Rebuild with: cmake -DCGEM_ENABLE_CALIBRATION=ON ..\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "For normal simulation runs, set calibration_mode=0 in INPUT/params.txt\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "================================================================================\n");
    fprintf(stderr, "\n");
    
    return -1;
}
