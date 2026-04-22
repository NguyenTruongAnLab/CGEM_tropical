/**
 * @file test_params_override_biogeo_rates.c
 * @brief Test that biogeochemical rate constants can be overridden via params.txt
 * 
 * This test verifies Phase 4 requirement: key biogeochemical parameters should be
 * configurable via INPUT/params.txt without recompiling.
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

    // Create a temporary params file with override values
    FILE *fp = fopen("test_params_override.txt", "w");
    if (!fp) {
        fprintf(stderr, "FAIL: Could not create temporary params file\n");
        return 1;
    }

    // Write test parameters with non-default values
    fprintf(fp, "# Test parameter overrides\n");
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
    
    // Phase 4 overrides: biogeochemical rate constants
    fprintf(fp, "Pb_Phy1 = 3.0e-4\n");           // Override max photosynthetic rate
    fprintf(fp, "Pb_Phy2 = 3.5e-4\n");
    fprintf(fp, "alpha_Phy1 = 1.0e-6\n");        // Override photosynthetic efficiency
    fprintf(fp, "kmortality_Phy1 = 2.0e-7\n");   // Override mortality rate
    fprintf(fp, "kox = 2.0e-4\n");               // Override aerobic degradation rate
    fprintf(fp, "knit = 2.0e-4\n");              // Override nitrification rate
    fprintf(fp, "pCO2atmo = 420.0\n");           // Override atmospheric CO2 (µatm)
    
    fclose(fp);

    // Read parameters (should set the override flags)
    read_parameters("test_params_override.txt");

    // Assign rate constants (should NOT overwrite params-set values)
    assignBiogeochemicalRateConstants();

    // Verify overrides took effect
    int failed = 0;

    if (fabs(Pb[Phy1] - 3.0e-4) > 1e-10) {
        fprintf(stderr, "FAIL: Pb[Phy1] not overridden (expected 3.0e-4, got %.10e)\n", Pb[Phy1]);
        failed = 1;
    }

    if (fabs(Pb[Phy2] - 3.5e-4) > 1e-10) {
        fprintf(stderr, "FAIL: Pb[Phy2] not overridden (expected 3.5e-4, got %.10e)\n", Pb[Phy2]);
        failed = 1;
    }

    if (fabs(alpha[Phy1] - 1.0e-6) > 1e-12) {
        fprintf(stderr, "FAIL: alpha[Phy1] not overridden (expected 1.0e-6, got %.10e)\n", alpha[Phy1]);
        failed = 1;
    }

    if (fabs(kmortality[Phy1] - 2.0e-7) > 1e-13) {
        fprintf(stderr, "FAIL: kmortality[Phy1] not overridden (expected 2.0e-7, got %.10e)\n", kmortality[Phy1]);
        failed = 1;
    }

    if (fabs(kox - 2.0e-4) > 1e-10) {
        fprintf(stderr, "FAIL: kox not overridden (expected 2.0e-4, got %.10e)\n", kox);
        failed = 1;
    }

    if (fabs(knit - 2.0e-4) > 1e-10) {
        fprintf(stderr, "FAIL: knit not overridden (expected 2.0e-4, got %.10e)\n", knit);
        failed = 1;
    }

    if (fabs(pCO2atmo - 420.0) > 0.1) {
        fprintf(stderr, "FAIL: pCO2atmo not overridden (expected 420.0, got %.2f)\n", pCO2atmo);
        failed = 1;
    }

    // Clean up
    remove("test_params_override.txt");

    if (failed) {
        return 1;
    }

    printf("PASS: All biogeochemical rate overrides applied correctly\n");
    return 0;
}
