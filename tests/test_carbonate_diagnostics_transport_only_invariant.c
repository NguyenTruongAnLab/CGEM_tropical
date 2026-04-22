/**
 * @file test_carbonate_diagnostics_transport_only_invariant.c
 * @brief Test that carbonate diagnostics (pH/CO2/pCO2) are computed in transport-only mode
 * 
 * This test verifies Phase 5 requirement: when enable_carbonate_diagnostics=1
 * and enable_biogeochemical_reactions=0, the model should compute pH/CO2/pCO2
 * from transported DIC+AT without modifying mass (mass-conservative diagnostic).
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "variables.h"

extern void Biogeo(int t);
extern void assignBiogeochemicalRateConstants(void);

int main(void)
{
    debug_level = 0;

    // 3-cell domain for better diagnosis
    M = 3;
    M1 = 2;
    M2 = 1;
    M3 = 0;

    DELXI = 1000;
    EL = 3000;
    WARMUP = 0;
    TS = 1;
    DELTI = 3600; // 1 hour timestep

    // CRITICAL: Transport-only mode with carbonate diagnostics enabled
    enable_biogeochemical_reactions = 0;
    enable_carbonate_diagnostics = 1;

    assignBiogeochemicalRateConstants();

    // Provide forcing data
    static double temp_time[2] = {0.0, 10.0};
    static double temp_val[2] = {28.0, 28.0};
    forcingData[FORCING_TEMPERATURE].time = temp_time;
    forcingData[FORCING_TEMPERATURE].data = temp_val;
    forcingData[FORCING_TEMPERATURE].dataSize = 2;

    static double light_time[2] = {0.0, 10.0};
    static double light_val[2] = {1000.0, 1000.0};
    forcingData[FORCING_LIGHT].time = light_time;
    forcingData[FORCING_LIGHT].data = light_val;
    forcingData[FORCING_LIGHT].dataSize = 2;

    static double wind_time[2] = {0.0, 10.0};
    static double wind_val[2] = {5.0, 5.0};
    forcingData[FORCING_WIND].time = wind_time;
    forcingData[FORCING_WIND].data = wind_val;
    forcingData[FORCING_WIND].dataSize = 2;

    // Initialize hydrodynamics
    for (int i = 1; i <= M; i++) {
        waterDepth[i] = 10.0;
        width[i] = 100.0;
        Chezy[i] = 50.0;
        velocity[i] = 0.5;
    }

    // Set up initial conditions with varying DIC/AT
    for (int i = 1; i <= M; i++) {
        v[Sal].c[i] = 5.0 + i * 2.0;  // Varying salinity
        v[Phy1].c[i] = 10.0;
        v[Phy2].c[i] = 5.0;
        v[Si].c[i] = 50.0;
        v[PO4].c[i] = 5.0;
        v[NO3].c[i] = 50.0;
        v[NH4].c[i] = 10.0;
        v[TOC].c[i] = 300.0;
        v[O2].c[i] = 160.0;
        v[DIC].c[i] = 1800.0 + i * 200.0;  // Varying DIC
        v[AT].c[i] = 2000.0 + i * 200.0;   // Varying AT
        v[SPM].c[i] = 10.0;
        
        // Initialize diagnostics to invalid values
        v[PH].c[i] = 0.0;
        v[CO2].c[i] = 0.0;
        v[pCO2].c[i] = 0.0;
    }

    // Store initial DIC/AT
    double DIC_initial[MAXM];
    double AT_initial[MAXM];
    for (int i = 1; i <= M; i++) {
        DIC_initial[i] = v[DIC].c[i];
        AT_initial[i] = v[AT].c[i];
    }

    // Run biogeochemistry (should only compute diagnostics, not modify mass)
    Biogeo(0);

    int failed = 0;

    // Verify that DIC and AT are unchanged (mass-conservative)
    for (int i = 1; i <= M; i++) {
        double DIC_change = fabs(v[DIC].c[i] - DIC_initial[i]);
        double AT_change = fabs(v[AT].c[i] - AT_initial[i]);
        
        if (DIC_change > 1e-6) {
            fprintf(stderr, "FAIL: DIC changed in transport-only mode at cell %d (Δ=%.10f)\n", 
                    i, DIC_change);
            failed = 1;
        }
        
        if (AT_change > 1e-6) {
            fprintf(stderr, "FAIL: AT changed in transport-only mode at cell %d (Δ=%.10f)\n", 
                    i, AT_change);
            failed = 1;
        }
    }

    // Verify that carbonate diagnostics were computed (not zero or NaN)
    for (int i = 1; i <= M; i++) {
        if (!isfinite(v[PH].c[i]) || v[PH].c[i] <= 0.0 || v[PH].c[i] > 14.0) {
            fprintf(stderr, "FAIL: Invalid pH at cell %d (pH=%.6f)\n", i, v[PH].c[i]);
            failed = 1;
        }
        
        if (!isfinite(v[CO2].c[i]) || v[CO2].c[i] < 0.0) {
            fprintf(stderr, "FAIL: Invalid CO2 at cell %d (CO2=%.6f)\n", i, v[CO2].c[i]);
            failed = 1;
        }
        
        if (!isfinite(v[pCO2].c[i]) || v[pCO2].c[i] < 0.0) {
            fprintf(stderr, "FAIL: Invalid pCO2 at cell %d (pCO2=%.6f)\n", i, v[pCO2].c[i]);
            failed = 1;
        }
        
        // Verify pH is in reasonable estuarine range
        if (v[PH].c[i] < 6.0 || v[PH].c[i] > 9.0) {
            fprintf(stderr, "FAIL: pH out of reasonable estuarine range at cell %d (pH=%.6f)\n", 
                    i, v[PH].c[i]);
            failed = 1;
        }
        
        printf("Cell %d: pH=%.3f, CO2=%.1f mmol/m³, pCO2=%.0f µatm (DIC=%.1f, AT=%.1f)\n",
               i, v[PH].c[i], v[CO2].c[i], v[pCO2].c[i], v[DIC].c[i], v[AT].c[i]);
    }

    if (failed) {
        return 1;
    }

    printf("PASS: Carbonate diagnostics computed correctly in transport-only mode without modifying DIC/AT\n");
    return 0;
}
