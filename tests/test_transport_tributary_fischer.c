#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "variables.h"
#include "transport.h"

static void init_minimal_grid(void)
{
    debug_level = 0;
    enable_flux_output = 0;
    enable_reaction_output = 0;
    calibration_mode = 0;

    WARMUP = 0;
    TS = 1;
    DELTI = 10;
    DELXI = 1000;

    M = 11;
    M1 = 10;
    M2 = 9;
    M3 = 11;

    for (int i = 1; i <= M; ++i) {
        totalArea[i] = 1000.0;
        width[i] = 100.0;
        disp[i] = 0.0;
        velocity[i] = 0.0;
    }
}

static int nearly_equal(double a, double b, double tol)
{
    return fabs(a - b) <= tol;
}

int main(void)
{
    init_minimal_grid();

    // Configure a single tributary in Fischer mode.
    tributaryEnabled = 1;
    numTributaries = 1;
    tributary_mass_injection_mode = TRIB_INJECTION_FISCHER;

    tributaries = (Tributary *)calloc((size_t)numTributaries, sizeof(Tributary));
    if (!tributaries) {
        fprintf(stderr, "FAIL: could not allocate tributaries\n");
        return 1;
    }

    // Place tributary at raw cell index 5 (odd continuity point).
    snprintf(tributaries[0].name, sizeof(tributaries[0].name), "%s", "TestTrib");
    tributaries[0].cellIndex = 5;

    // Provide discharge time series (m3/s). Must be positive for tributaries.
    double discharge_time[2] = {0.0, 1.0};
    double discharge_val[2] = {10.0, 10.0};
    tributaries[0].dischargeTime = discharge_time;
    tributaries[0].discharge = discharge_val;
    tributaries[0].dischargeDataSize = 2;

    // Provide tracer concentration time series at the tributary (mmol/m3).
    const int s = NO3;
    double conc_time[2] = {0.0, 1.0};
    double conc_val[2] = {100.0, 100.0};
    tributaries[0].chemicalData[s].timeArray = conc_time;
    tributaries[0].chemicalData[s].dataArray = conc_val;
    tributaries[0].chemicalData[s].dataSize = 2;

    // Initial condition: zero tracer everywhere.
    double c[MAXM + 1];
    for (int i = 1; i <= M; ++i) c[i] = 0.0;

    // Expected single-step Fischer relaxation on the continuity control volume.
    // For interior odd index, applyTributarySourceTerms uses dx_control = 2*DELXI.
    const int cell = 5;
    const double A_cell = totalArea[cell];
    const double dx_control = 2.0 * (double)DELXI;
    const double V = A_cell * dx_control;

    const double Q_trib = discharge_val[0];
    const double C_trib = conc_val[0];
    const double dt = (double)DELTI;

    // From code comments: C_new = C_old + alpha*(C_trib - C_old), alpha = 1 - exp(-(Q/V)*dt)
    const double alpha = 1.0 - exp(-(Q_trib / V) * dt);
    const double expected = 0.0 + alpha * (C_trib - 0.0);

    applyTributarySourceTerms(c, s, 0);

    if (!nearly_equal(c[cell], expected, 1e-10)) {
        fprintf(stderr,
                "FAIL: Fischer tributary relaxation mismatch at cell %d (got %.15g expected %.15g; alpha=%.15g)\n",
                cell, c[cell], expected, alpha);
        free(tributaries);
        tributaries = NULL;
        return 1;
    }

    // Sanity: no other cells should be modified for Fischer single-cell injection.
    for (int i = 1; i <= M; ++i) {
        if (i == cell) continue;
        if (!nearly_equal(c[i], 0.0, 1e-14)) {
            fprintf(stderr, "FAIL: unexpected concentration change at i=%d (%.15g)\n", i, c[i]);
            free(tributaries);
            tributaries = NULL;
            return 1;
        }
    }

    free(tributaries);
    tributaries = NULL;

    return 0;
}
