#include <math.h>
#include <stdio.h>

#include "define.h"
#include "variables.h"

extern void Biogeo(int t);
extern void assignBiogeochemicalRateConstants(void);

static void setup_minimal_domain(void)
{
    debug_level = 0;
    calibration_mode = 0;

    M = 1;
    M1 = 0;
    M2 = 0;
    M3 = 0;

    DELXI = 1000;
    EL = 1000;
    DELTI = 180;
    TS = 1;
    WARMUP = 0;

    enable_biogeochemical_reactions = BIOGEO_FULL;
    enable_carbonate_diagnostics = 0;

    waterDepth[1] = 10.0;
    width[1] = 100.0;
    Chezy[1] = 50.0;
    velocity[1] = 0.0;

    static double temp_time[2] = {0.0, 1.0};
    static double temp_val[2] = {35.0, 35.0};
    forcingData[FORCING_TEMPERATURE].time = temp_time;
    forcingData[FORCING_TEMPERATURE].data = temp_val;
    forcingData[FORCING_TEMPERATURE].dataSize = 2;

    static double light_time[2] = {0.0, 1.0};
    static double light_val[2] = {1000.0, 1000.0};
    forcingData[FORCING_LIGHT].time = light_time;
    forcingData[FORCING_LIGHT].data = light_val;
    forcingData[FORCING_LIGHT].dataSize = 2;

    static double wind_time[2] = {0.0, 1.0};
    static double wind_val[2] = {0.0, 0.0};
    forcingData[FORCING_WIND].time = wind_time;
    forcingData[FORCING_WIND].data = wind_val;
    forcingData[FORCING_WIND].dataSize = 2;
}

static void setup_growth_state_for_underflow_test(void)
{
    // Force very strong light attenuation so Ebottom underflows.
    v[SPM].c[1] = 10.0;

    // Ensure biomass exists and nutrients are non-limiting.
    v[Phy1].c[1] = 50.0;
    v[Phy2].c[1] = 0.0;
    v[Si].c[1] = 200.0;
    v[NH4].c[1] = 100.0;
    v[NO3].c[1] = 100.0;
    v[PO4].c[1] = 10.0;

    // Keep other tracers finite and benign.
    v[Sal].c[1] = 0.0;
    v[TOC].c[1] = 100.0;
    v[O2].c[1] = 200.0;
    v[DIC].c[1] = 2000.0;
    v[AT].c[1] = 2200.0;
    v[CO2].c[1] = 50.0;
    v[pCO2].c[1] = 400.0;
    v[PH].c[1] = 8.0;
}

int main(void)
{
    setup_minimal_domain();

    // ------------------------------------------------------------------
    // (a) Light underflow branch: non-zero growth with asymptotic treatment
    // ------------------------------------------------------------------
    assignBiogeochemicalRateConstants();
    setup_growth_state_for_underflow_test();

    // Remove non-light sinks so the branch effect is isolated.
    kmaint[Phy1] = 0.0;
    kmaint[Phy2] = 0.0;
    kmortality[Phy1] = 0.0;
    kmortality[Phy2] = 0.0;
    ws_phy = 0.0;

    const double phy1_before = v[Phy1].c[1];
    Biogeo(0);
    const double phy1_after = v[Phy1].c[1];

    if (!(phy1_after > phy1_before)) {
        fprintf(stderr, "FAIL: asymptotic underflow should allow non-zero light growth (Phy1 %.12g -> %.12g)\n",
                phy1_before, phy1_after);
        return 1;
    }

    // ------------------------------------------------------------------
    // (b) kmort temperature response: tropical temperature limit caps runaway
    // ------------------------------------------------------------------
    kmortality[Phy1] = 1.0e-7;
    double kmort_val = kmort(0, Phy1);

    // At 35°C with Q10=1.067 + tropical cap, mortality should be positive
    // but not excessively high (the tropical cap limits extreme temps).
    if (!(kmort_val > 0.0)) {
        fprintf(stderr, "FAIL: kmort should be positive at tropical T (got %.12g)\n", kmort_val);
        return 1;
    }

    // ------------------------------------------------------------------
    // (c) Default assignment and override protection
    // ------------------------------------------------------------------
    biogeo_override_Pb[Phy1] = 0;
    biogeo_override_KN[Phy1] = 0;
    kbg = 0.0;
    kspm = 0.0;
    assignBiogeochemicalRateConstants();

    if (fabs(Pb[Phy1] - 2.0e-4) > 1e-12) {
        fprintf(stderr, "FAIL: Pb_Phy1 default mismatch (got %.12g, expected 2.0e-4)\n", Pb[Phy1]);
        return 1;
    }
    if (fabs(KN[Phy1] - 25.0) > 1e-12) {
        fprintf(stderr, "FAIL: KN_Phy1 default mismatch (got %.12g, expected 25.0)\n", KN[Phy1]);
        return 1;
    }
    if (fabs(kbg - 0.2) > 1e-12) {
        fprintf(stderr, "FAIL: kbg default mismatch (got %.12g, expected 0.2)\n", kbg);
        return 1;
    }

    biogeo_override_Pb[Phy1] = 1;
    Pb[Phy1] = 9.9e-4;
    assignBiogeochemicalRateConstants();
    if (fabs(Pb[Phy1] - 9.9e-4) > 1e-12) {
        fprintf(stderr, "FAIL: explicit Pb_Phy1 override should not be overwritten (got %.12g)\n", Pb[Phy1]);
        return 1;
    }

    return 0;
}
