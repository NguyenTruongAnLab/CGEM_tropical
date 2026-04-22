#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "variables.h"

// biogeo.c
extern void Biogeo(int t);
// init.c
extern void assignBiogeochemicalRateConstants(void);

static int isfinite_pos(double x) {
    return isfinite(x) && x > 0.0;
}

int main(void)
{
    debug_level = 0;

    // Minimal 1-cell domain.
    DELTI = 10;
    TS = 1;
    WARMUP = 0;

    M = 1;
    M1 = 0;
    M2 = 0;
    M3 = 0;

    // Ensure carbonate diagnostics run without reactions.
    enable_biogeochemical_reactions = 0;
    enable_carbonate_diagnostics = 1;

    // Provide reasonable defaults for constants used by carbonate diagnostics.
    assignBiogeochemicalRateConstants();

    waterDepth[1] = 10.0;
    width[1] = 100.0;
    velocity[1] = 0.0;

    // Salinity controls equilibrium constants.
    v[Sal].c[1] = 10.0;

    // Choose a typical estuarine carbonate state.
    v[DIC].c[1] = 2000.0; // mmol/m3
    v[AT].c[1]  = 2200.0; // mmol/m3

    // Baseline nutrients.
    v[NH4].c[1] = 10.0;
    v[NO3].c[1] = 0.0;

    Biogeo(0);
    const double pH_base = v[PH].c[1];
    if (!isfinite_pos(pH_base)) {
        fprintf(stderr, "FAIL: baseline pH not finite/positive (%g)\n", pH_base);
        return 1;
    }

    // Increase NO3 at fixed DIC/AT. With nutrient-alkalinity consistency (Phase 5),
    // higher NO3 reduces TA available to carbonate in the balance (TA includes -NO3),
    // so the carbonate system must become more basic (higher pH) to match the same total TA.
    v[NO3].c[1] = 200.0;
    Biogeo(0);
    const double pH_high_no3 = v[PH].c[1];
    if (!isfinite_pos(pH_high_no3)) {
        fprintf(stderr, "FAIL: high-NO3 pH not finite/positive (%g)\n", pH_high_no3);
        return 1;
    }

    if (!(pH_high_no3 > pH_base)) {
        fprintf(stderr, "FAIL: expected pH to increase when NO3 increases at fixed DIC/AT (base=%g highNO3=%g)\n", pH_base, pH_high_no3);
        return 1;
    }

    // Increase NH4 at fixed DIC/AT (and low NO3). With Phase 5, higher NH4 increases the non-carbonate
    // alkalinity term (+NH4), so carbonate alkalinity can be lower and pH decreases.
    v[NO3].c[1] = 0.0;
    v[NH4].c[1] = 200.0;
    Biogeo(0);
    const double pH_high_nh4 = v[PH].c[1];
    if (!isfinite_pos(pH_high_nh4)) {
        fprintf(stderr, "FAIL: high-NH4 pH not finite/positive (%g)\n", pH_high_nh4);
        return 1;
    }

    if (!(pH_high_nh4 < pH_base)) {
        fprintf(stderr, "FAIL: expected pH to decrease when NH4 increases at fixed DIC/AT (base=%g highNH4=%g)\n", pH_base, pH_high_nh4);
        return 1;
    }

    return 0;
}
