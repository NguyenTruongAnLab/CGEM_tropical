#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "variables.h"
#include "transport.h"

static double compute_total_mass_odd_centers(const double *c)
{
    double total = 0.0;

    for (int i = 1; i <= M; i += 2) {
        const double A = totalArea[i];
        if (!(A > 0.0) || !isfinite(A)) {
            continue;
        }

        // Odd-index control volumes span ~2*DELXI, except half-cells at the boundaries.
        double dx = 2.0 * (double)DELXI;
        if (i == 1 || i == M) {
            dx = (double)DELXI;
        }

        const double Ci = (isfinite(c[i]) && c[i] > 0.0) ? c[i] : 0.0;
        total += A * Ci * dx;
    }

    return total;
}

int main(void)
{
    // Keep the test deterministic and independent from INPUT/ forcing.
    debug_level = 0;
    enable_flux_output = 0;
    enable_reaction_output = 0;
    calibration_mode = 0;

    WARMUP = 0;
    TS = 1;
    DELTI = 10;
    DELXI = 1000;

    // Staggered grid convention used by transport.c TVD(): M is odd so M1 is even (face index).
    M = 11;
    M1 = 10;
    M2 = 9;
    M3 = 11;

    // Geometry and flow field
    for (int i = 1; i <= M; ++i) {
        totalArea[i] = 1000.0;
        width[i] = 100.0;
        disp[i] = 0.0;
        velocity[i] = 0.0;
    }

    // Closed domain for this unit test: zero flux at the boundary faces.
    // Use a uniform interior face velocity to advect mass without creating/removing it.
    for (int i = 2; i <= M1; i += 2) {
        velocity[i] = 0.2;
    }
    velocity[2] = 0.0;
    velocity[M1] = 0.0;

    double c[MAXM + 1];
    for (int i = 1; i <= M; ++i) {
        c[i] = 0.0;
    }

    // Non-uniform interior tracer (odd centers only)
    c[5] = 1.0;
    c[7] = 2.0;
    c[9] = 1.0;

    const double mass0 = compute_total_mass_odd_centers(c);

    const int steps = 200;
    for (int n = 0; n < steps; ++n) {
        TVD(c, Sal);
    }

    const double mass1 = compute_total_mass_odd_centers(c);

    const double denom = (mass0 > 0.0) ? mass0 : 1.0;
    const double rel_err = fabs(mass1 - mass0) / denom;

    // Double-precision roundoff should keep this extremely small.
    if (!(rel_err < 1e-10)) {
        fprintf(stderr, "FAIL: mass not conserved (mass0=%.15e mass1=%.15e rel=%.3e)\n", mass0, mass1, rel_err);
        return 1;
    }

    for (int i = 1; i <= M; i += 2) {
        if (!isfinite(c[i])) {
            fprintf(stderr, "FAIL: non-finite concentration at i=%d\n", i);
            return 1;
        }
        if (c[i] < -1e-14) {
            fprintf(stderr, "FAIL: negative concentration at i=%d (%g)\n", i, c[i]);
            return 1;
        }
    }

    return 0;
}
