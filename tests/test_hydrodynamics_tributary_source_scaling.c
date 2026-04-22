#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "variables.h"

// hydrodynamics.c does not expose a header; this function is non-static.
extern void Coeffa(int t);

static int nearly_equal(double a, double b, double tol)
{
    return fabs(a - b) <= tol;
}

int main(void)
{
    // Deterministic, independent from INPUT.
    debug_level = 0;
    calibration_mode = 0;

    // Minimal grid consistent with runtime convention: M computed even.
    DELTI = 10;
    DELXI = 1000;
    M = 10;
    M1 = M - 1;
    M2 = M - 2;
    M3 = M - 3;

    // Ensure forcing discharge is available for Discharge_ups() during Coeffa().
    static double q_time[2] = {0.0, 86400.0};
    static double q_val[2] = {100.0, 100.0};
    forcingData[FORCING_DISCHARGE].time = q_time;
    forcingData[FORCING_DISCHARGE].data = q_val;
    forcingData[FORCING_DISCHARGE].dataSize = 2;

    // Minimal geometry/hydro fields so Coeffa() doesn't encounter divisions by zero.
    for (int i = 1; i <= M; ++i) {
        width[i] = 100.0;
        Chezy[i] = 50.0;
        waterDepth[i] = 10.0;

        baseArea[i] = width[i] * waterDepth[i];
        freeArea[i] = 0.0;
        tempFreeArea[i] = freeArea[i];
        totalArea[i] = baseArea[i] + freeArea[i];

        velocity[i] = 0.0;
        tempVelocity[i] = 0.0;

        // Clear any prior coefficients.
        Z[i] = 0.0;
        for (int k = 0; k < 5; ++k) {
            C[i][k] = 0.0;
        }
    }

    // Configure a single tributary: Q_trib should appear as a continuity source term.
    tributaryEnabled = 1;
    numTributaries = 1;

    tributaries = (Tributary *)calloc((size_t)numTributaries, sizeof(Tributary));
    if (!tributaries) {
        fprintf(stderr, "FAIL: could not allocate tributaries\n");
        return 1;
    }

    // Put tributary at cell 5 (odd, continuity point).
    tributaries[0].cellIndex = 5;
    snprintf(tributaries[0].name, sizeof(tributaries[0].name), "%s", "TestTrib");

    static double trib_time[2] = {0.0, 86400.0};
    static double trib_q[2] = {10.0, 10.0};
    tributaries[0].dischargeTime = trib_time;
    tributaries[0].discharge = trib_q;
    tributaries[0].dischargeDataSize = 2;

    // Run coefficient assembly at t=0.
    Coeffa(0);

    // Continuity rows are on odd indices in Coeffa() interior loop.
    const int j = 5;

    // With freeArea[j] set to 0, Z[j] should equal the lateral inflow source term [m^2/s].
    // Discrete continuity uses inv_2dx = 1/(2*DELXI) for flux divergence, so a tributary discharge
    // applied over the same control length should contribute Q_trib/(2*DELXI).
    const double expected = trib_q[0] / (2.0 * (double)DELXI);
    const double got = Z[j];

    if (!nearly_equal(got, expected, 1e-12)) {
        fprintf(stderr, "FAIL: tributary continuity source scaling mismatch at j=%d (got %.15g expected %.15g)\n", j, got, expected);
        free(tributaries);
        tributaries = NULL;
        return 1;
    }

    free(tributaries);
    tributaries = NULL;
    return 0;
}
