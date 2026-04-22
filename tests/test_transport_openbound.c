#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "variables.h"
#include "transport.h"

static void init_minimal_grid(void)
{
    // Deterministic, independent of INPUT.
    debug_level = 0;
    enable_flux_output = 0;
    enable_reaction_output = 0;
    calibration_mode = 0;

    WARMUP = 0;
    TS = 1;
    DELTI = 10;
    DELXI = 1000;

    // Odd M so M1 is even face index.
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

    // Boundary faces: mouth at index 2, upstream face at M1.
    // These are set per-test.
}

static int nearly_equal(double a, double b, double tol)
{
    return fabs(a - b) <= tol;
}

int main(void)
{
    init_minimal_grid();

    // Generic tracer boundary test (NO3): CFL-bounded upwind blending at mouth.
    {
        const int s = NO3;
        v[s].clb = 20.0; // mouth reservoir
        v[s].cub = 3.0;  // upstream reservoir (not under test)

        double c[MAXM + 1];
        for (int i = 1; i <= M; ++i) c[i] = 5.0;
        c[3] = 6.0;

        // FLOOD (inflow at mouth): source is clb, blended by alpha = u*dt/dx.
        velocity[2] = 0.1;
        disp[2] = 0.0; // isolate advection enforcement

        Openbound(c, s, 0);
        {
            const double alpha = fmin(velocity[2] * (double)DELTI / (double)DELXI, 1.0);
            const double expected = 5.0 + alpha * (v[s].clb - 5.0);
            if (!nearly_equal(c[1], expected, 1e-12)) {
                fprintf(stderr, "FAIL: expected generic tracer flood BC CFL blend (got %.15g expected %.15g)\n", c[1], expected);
                return 1;
            }
        }

        // EBB (outflow at mouth): source is interior cell 3, blended by alpha = |u|*dt/dx.
        for (int i = 1; i <= M; ++i) c[i] = 7.0;
        c[3] = 9.0;
        velocity[2] = -0.1;
        disp[2] = 0.0; // isolate advection

        Openbound(c, s, 0);
        {
            const double alpha = fmin(-velocity[2] * (double)DELTI / (double)DELXI, 1.0);
            const double expected = 7.0 + alpha * (9.0 - 7.0);
            if (!nearly_equal(c[1], expected, 1e-12)) {
                fprintf(stderr, "FAIL: expected generic tracer ebb BC CFL blend (got %.15g expected %.15g)\n", c[1], expected);
                return 1;
            }
        }

        if (!isfinite(c[1]) || c[1] < 0.0) {
            fprintf(stderr, "FAIL: invalid mouth concentration after ebb blend (%.15g)\n", c[1]);
            return 1;
        }
    }

    // Salinity expected behavior (Phase 3 target): during ebb (outflow), do not pin to clb.
    {
        const int s = Sal;
        v[s].clb = 30.0;
        v[s].cub = 0.0;

        double c[MAXM + 1];
        for (int i = 1; i <= M; ++i) c[i] = 0.0;

        // EBB (outflow): velocity negative.
        velocity[2] = -0.1;
        disp[2] = 0.0;

        Openbound(c, s, 0);

        // With zero initial salinity and no dispersive exchange, ebb should not introduce ocean salinity.
        if (!nearly_equal(c[1], 0.0, 1e-12)) {
            fprintf(stderr, "FAIL: expected salinity ebb BC to not pin to clb (got %.15g)\n", c[1]);
            return 1;
        }
    }

    // Salinity flood follows the same CFL-upwind blending (parameter-free).
    {
        const int s = Sal;
        v[s].clb = 30.0;
        v[s].cub = 0.0;

        double c[MAXM + 1];
        for (int i = 1; i <= M; ++i) c[i] = 5.0;
        c[3] = 6.0;

        velocity[2] = 0.5; // flood
        disp[2] = 0.0;

        Openbound(c, s, 0);

        const double alpha = fmin(velocity[2] * (double)DELTI / (double)DELXI, 1.0);
        const double expected = 5.0 + alpha * (v[s].clb - 5.0);
        if (!nearly_equal(c[1], expected, 1e-12)) {
            fprintf(stderr, "FAIL: expected salinity flood CFL blend (got %.15g expected %.15g)\n", c[1], expected);
            return 1;
        }
    }

    // Slack-tide continuity check: no artificial jump around u≈0 for salinity BC.
    {
        const int s = Sal;
        v[s].clb = 35.0;
        v[s].cub = 0.0;

        double c0[MAXM + 1], cplus[MAXM + 1], cminus[MAXM + 1];
        for (int i = 1; i <= M; ++i) {
            c0[i] = 10.0;
            cplus[i] = 10.0;
            cminus[i] = 10.0;
        }
        c0[3] = cplus[3] = cminus[3] = 12.0;

        velocity[2] = 0.0;
        Openbound(c0, s, 0);
        if (!nearly_equal(c0[1], 10.0, 1e-12)) {
            fprintf(stderr, "FAIL: expected zero mouth velocity to keep salinity unchanged (got %.15g)\n", c0[1]);
            return 1;
        }

        velocity[2] = 1e-6;
        Openbound(cplus, s, 0);
        velocity[2] = -1e-6;
        Openbound(cminus, s, 0);

        if (fabs(cplus[1] - 10.0) > 1e-6 || fabs(cminus[1] - 10.0) > 1e-6) {
            fprintf(stderr, "FAIL: expected near-slack salinity BC to remain near prior state (u=±1e-6, got + %.15g, - %.15g)\n", cplus[1], cminus[1]);
            return 1;
        }
        if (fabs(cplus[1] - cminus[1]) > 1e-6) {
            fprintf(stderr, "FAIL: expected no artificial salinity jump across slack tide (plus %.15g minus %.15g)\n", cplus[1], cminus[1]);
            return 1;
        }
    }

    return 0;
}
