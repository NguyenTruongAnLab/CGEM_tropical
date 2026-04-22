/**
 * @file test_missing_param_behavior_is_documented.c
 * @brief Test that when params are missing, compiled defaults are used
 * 
 * This test verifies Phase 4 requirement: defaults must match current
 * compiled behavior when keys are absent from params.txt.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "variables.h"
#include "file.h"

extern void assignBiogeochemicalRateConstants(void);

int main(void)
{
    debug_level = 0;
    M = 1;
    DELTI = 180;
    WARMUP = 0;

    // Create a minimal params file WITHOUT biogeochemical overrides
    FILE *fp = fopen("test_params_minimal.txt", "w");
    if (!fp) {
        fprintf(stderr, "FAIL: Could not create temporary params file\n");
        return 1;
    }

    // Write minimal required parameters only
    fprintf(fp, "# Minimal test parameters\n");
    fprintf(fp, "MAXT = 1\n");
    fprintf(fp, "WARMUP = 0\n");
    fprintf(fp, "DELTI = 180\n");
    fprintf(fp, "TS = 1\n");
    fprintf(fp, "DELXI = 1000\n");
    fprintf(fp, "EL = 1000\n");
    fprintf(fp, "AMPL = 1.0\n");
    fprintf(fp, "num_segments = 1\n");
    fprintf(fp, "index_1 = 0\n");
    fprintf(fp, "B1 = 1000.0\n");
    fprintf(fp, "LC1 = 50000.0\n");
    fprintf(fp, "Chezy1 = 50.0\n");
    fprintf(fp, "Rs1 = 1.0\n");
    
    fclose(fp);

    // Read parameters (no biogeochemical overrides present)
    read_parameters("test_params_minimal.txt");

    // Assign rate constants (should use compiled defaults)
    assignBiogeochemicalRateConstants();

    // Verify compiled defaults are used (from init.c::assignBiogeochemicalRateConstants)
    int failed = 0;

    // Check documented defaults match compiled values
    const double expected_Pb_Phy1 = 2.0e-4;
    const double expected_Pb_Phy2 = 2.6e-4;
    const double expected_alpha_Phy1 = 8.0e-7;
    const double expected_kmortality_Phy1 = 1.4e-7;
    const double expected_kox = 1.5e-4;
    const double expected_knit = 1.5e-4;
    const double expected_pCO2atmo = 415.0;  // Updated to 2024 global mean

    if (fabs(Pb[Phy1] - expected_Pb_Phy1) > 1e-10) {
        fprintf(stderr, "FAIL: Default Pb[Phy1] mismatch (expected %.10e, got %.10e)\n", 
                expected_Pb_Phy1, Pb[Phy1]);
        failed = 1;
    }

    if (fabs(Pb[Phy2] - expected_Pb_Phy2) > 1e-10) {
        fprintf(stderr, "FAIL: Default Pb[Phy2] mismatch (expected %.10e, got %.10e)\n", 
                expected_Pb_Phy2, Pb[Phy2]);
        failed = 1;
    }

    if (fabs(alpha[Phy1] - expected_alpha_Phy1) > 1e-12) {
        fprintf(stderr, "FAIL: Default alpha[Phy1] mismatch (expected %.10e, got %.10e)\n", 
                expected_alpha_Phy1, alpha[Phy1]);
        failed = 1;
    }

    if (fabs(kmortality[Phy1] - expected_kmortality_Phy1) > 1e-13) {
        fprintf(stderr, "FAIL: Default kmortality[Phy1] mismatch (expected %.10e, got %.10e)\n", 
                expected_kmortality_Phy1, kmortality[Phy1]);
        failed = 1;
    }

    if (fabs(kox - expected_kox) > 1e-10) {
        fprintf(stderr, "FAIL: Default kox mismatch (expected %.10e, got %.10e)\n", 
                expected_kox, kox);
        failed = 1;
    }

    if (fabs(knit - expected_knit) > 1e-10) {
        fprintf(stderr, "FAIL: Default knit mismatch (expected %.10e, got %.10e)\n", 
                expected_knit, knit);
        failed = 1;
    }

    if (fabs(pCO2atmo - expected_pCO2atmo) > 0.1) {
        fprintf(stderr, "FAIL: Default pCO2atmo mismatch (expected %.2f, got %.2f)\n", 
                expected_pCO2atmo, pCO2atmo);
        failed = 1;
    }

    // Clean up
    remove("test_params_minimal.txt");

    if (failed) {
        return 1;
    }

    printf("PASS: All default biogeochemical rates match compiled behavior\n");
    return 0;
}
