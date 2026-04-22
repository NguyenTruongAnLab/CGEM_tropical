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

    // 1-cell domain to isolate kinetics.
    M = 1;
    M1 = 0;
    M2 = 0;
    M3 = 0;

    // Needed by windspeed() spatial decay.
    DELXI = 1000;
    EL = 1000;

    WARMUP = 0;
    TS = 1;

    // Use a deliberately large timestep to trigger overshoot if rates are not limited.
    DELTI = 6 * 3600; // 6 hours

    enable_biogeochemical_reactions = 1;
    enable_carbonate_diagnostics = 0; // not needed for this test

    // Ensure biogeo constants are set.
    assignBiogeochemicalRateConstants();

    // Provide forcing data required by waterT(t) and I0(t) in biogeout.c.
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

    // Minimal hydro/geometry.
    waterDepth[1] = 10.0;
    width[1] = 100.0;
    Chezy[1] = 50.0;
    velocity[1] = 0.0;

    // Initial concentrations chosen to be small enough that naive explicit updates will go negative.
    v[Sal].c[1] = 0.0;

    v[Phy1].c[1] = 50.0;
    v[Phy2].c[1] = 50.0;

    v[Si].c[1] = 1.0;
    v[PO4].c[1] = 0.2;
    v[NO3].c[1] = 1.0;
    v[NH4].c[1] = 1.0;

    v[TOC].c[1] = 5.0;
    v[O2].c[1] = 2.0;

    v[DIC].c[1] = 2000.0;
    v[AT].c[1] = 2200.0;

    Biogeo(0);

    // Also assert non-negativity as a sanity check.
    for (int s = 0; s < CHEM_COUNT; ++s) {
        if (!isfinite(v[s].c[1]) || v[s].c[1] < 0.0) {
            fprintf(stderr, "FAIL: invalid concentration after Biogeo for species %d (%g)\n", s, v[s].c[1]);
            return 1;
        }
    }

    return 0;
}
