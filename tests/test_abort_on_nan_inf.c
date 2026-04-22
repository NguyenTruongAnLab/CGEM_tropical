/**
 * @file test_abort_on_nan_inf.c
 * @brief Test that the model aborts when any tracer becomes NaN or Inf.
 * 
 * This is a "death test": we expect the program to abort with a fatal error message,
 * not to continue with invalid floating-point values.
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

    M = 1;
    M1 = 0;
    M2 = 0;
    M3 = 0;

    DELXI = 1000;
    EL = 1000;

    WARMUP = 0;
    TS = 1;
    DELTI = 3600;

    enable_biogeochemical_reactions = 1;
    enable_carbonate_diagnostics = 0;

    assignBiogeochemicalRateConstants();

    static double temp_time[2] = {0.0, 10.0};
    static double temp_val[2] = {28.0, 28.0};
    forcingData[FORCING_TEMPERATURE].time = temp_time;
    forcingData[FORCING_TEMPERATURE].data = temp_val;
    forcingData[FORCING_TEMPERATURE].dataSize = 2;

    static double light_time[2] = {0.0, 10.0};
    static double light_val[2] = {3000.0, 3000.0};
    forcingData[FORCING_LIGHT].time = light_time;
    forcingData[FORCING_LIGHT].data = light_val;
    forcingData[FORCING_LIGHT].dataSize = 2;

    static double wind_time[2] = {0.0, 10.0};
    static double wind_val[2] = {0.0, 0.0};
    forcingData[FORCING_WIND].time = wind_time;
    forcingData[FORCING_WIND].data = wind_val;
    forcingData[FORCING_WIND].dataSize = 2;

    waterDepth[1] = 10.0;
    width[1] = 100.0;
    Chezy[1] = 50.0;
    velocity[1] = 0.0;

    // Set up normal initial conditions
    v[Sal].c[1] = 0.0;
    v[Phy1].c[1] = 50.0;
    v[Phy2].c[1] = 50.0;
    v[Si].c[1] = 30.0;
    v[PO4].c[1] = 1.0;
    v[NO3].c[1] = 20.0;
    v[NH4].c[1] = 5.0;
    v[TOC].c[1] = 300.0;
    v[O2].c[1] = 160.0;
    v[DIC].c[1] = 2000.0;
    v[AT].c[1] = 2200.0;

    // Inject NaN into a tracer - should be detected and cause abort
    v[NO3].c[1] = NAN;

    // This should detect NaN and abort.
    Biogeo(0);

    // If we get here, the test FAILED (we expected an abort).
    fprintf(stderr, "FAIL: Expected abort on NaN concentration, but program continued.\n");
    return 1;
}
