/**
 * @file test_transport_boundary_cv_length.c
 * @brief Unit test: Verify boundary control-volume length convention (DELXI)
 * 
 * Phase 2 TDD Requirement:
 * Pin down the staggered-grid boundary CV length assumption for consistency:
 * - Interior odd CV length = 2*DELXI
 * - Boundary half-cell length = DELXI (not 0.5*DELXI)
 * 
 * This test validates:
 * 1. Openbound() outflow removal rate (pure advection scenario)
 * 2. Fischer tributary injection at boundary cell
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

    M = 7;
    M1 = 6;
    M2 = 5;
    M3 = 7;

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

/**
 * @brief Test 1: Openbound() outflow removal rate assumes boundary CV length = DELXI
 * 
 * Scenario: Pure advection outflow at mouth (Q_mouth < 0, D = 0).
 * After boundary update, remaining mass fraction should be:
 *   remaining_fraction = 1 - (|Q| * dt) / (A * DELXI)
 * 
 * This fails if the code uses volume = A * DELXI * 0.5 (current state).
 */
static int test_openbound_outflow_removal_rate(void)
{
    init_minimal_grid();

    // Configure pure outflow at mouth (no dispersion)
    const int mouth_face = 2;
    const double A_mouth = totalArea[1];
    const double Q_out = 100.0;  // m³/s, outflow magnitude
    // In Openbound(), mouth inflow is Q_mouth > 0; outflow is Q_mouth < 0.
    velocity[mouth_face] = -Q_out / totalArea[mouth_face];
    disp[mouth_face] = 0.0;

    // Initial concentration at mouth cell
    double c[MAXM + 1];
    for (int i = 1; i <= M; ++i) c[i] = 0.0;
    c[1] = 10.0;  // mmol/m³

    // Boundary condition (not used during outflow)
    v[Sal].clb = 0.0;
    v[Sal].cub = 0.0;

    // Expected calculation using boundary CV length = DELXI (Option A):
    const double volume_mouth = A_mouth * (double)DELXI;  // NOT 0.5*DELXI
    const double dt = (double)DELTI;
    const double outflow_volume = Q_out * dt;
    const double remove_fraction = fmin(outflow_volume / volume_mouth, 1.0);
    const double expected = c[1] * (1.0 - remove_fraction);

    // Run Openbound
    Openbound(c, Sal, 0);

    if (!nearly_equal(c[1], expected, 1e-10)) {
        fprintf(stderr,
                "FAIL: Openbound outflow removal rate mismatch (got %.15g expected %.15g; remove_fraction=%.15g)\n",
                c[1], expected, remove_fraction);
        fprintf(stderr,
                "  → This indicates boundary volume uses 0.5*DELXI instead of DELXI\n");
        return 1;
    }

    return 0;
}

/**
 * @brief Test 2: Fischer tributary injection at boundary cell assumes dx_control = DELXI
 * 
 * Scenario: Tributary injecting into boundary cell (i=1 or i=M).
 * The relaxation mixing rate should use volume = A * DELXI (not A * 0.5*DELXI).
 */
static int test_fischer_tributary_boundary_cell(void)
{
    init_minimal_grid();

    tributaryEnabled = 1;
    numTributaries = 1;
    tributary_mass_injection_mode = TRIB_INJECTION_FISCHER;

    tributaries = (Tributary *)calloc(1, sizeof(Tributary));
    if (!tributaries) {
        fprintf(stderr, "FAIL: could not allocate tributaries\n");
        return 1;
    }

    // Place tributary at boundary cell (i=1, mouth)
    snprintf(tributaries[0].name, sizeof(tributaries[0].name), "%s", "BoundaryTrib");
    tributaries[0].cellIndex = 1;

    double discharge_time[2] = {0.0, 1.0};
    double discharge_val[2] = {5.0, 5.0};
    tributaries[0].dischargeTime = discharge_time;
    tributaries[0].discharge = discharge_val;
    tributaries[0].dischargeDataSize = 2;

    const int s = NO3;
    double conc_time[2] = {0.0, 1.0};
    double conc_val[2] = {50.0, 50.0};
    tributaries[0].chemicalData[s].timeArray = conc_time;
    tributaries[0].chemicalData[s].dataArray = conc_val;
    tributaries[0].chemicalData[s].dataSize = 2;

    double c[MAXM + 1];
    for (int i = 1; i <= M; ++i) c[i] = 0.0;

    const int cell = 1;
    const double A_cell = totalArea[cell];
    // Expected: boundary dx_control = DELXI (not 0.5*DELXI)
    const double dx_control = (double)DELXI;
    const double V = A_cell * dx_control;

    const double Q_trib = discharge_val[0];
    const double C_trib = conc_val[0];
    const double dt = (double)DELTI;

    const double alpha = 1.0 - exp(-(Q_trib / V) * dt);
    const double expected = 0.0 + alpha * (C_trib - 0.0);

    applyTributarySourceTerms(c, s, 0);

    if (!nearly_equal(c[cell], expected, 1e-10)) {
        fprintf(stderr,
                "FAIL: Fischer tributary at boundary cell mismatch (got %.15g expected %.15g; alpha=%.15g)\n",
                c[cell], expected, alpha);
        fprintf(stderr,
                "  → This indicates boundary dx_control uses 0.5*DELXI instead of DELXI\n");
        free(tributaries);
        tributaries = NULL;
        return 1;
    }

    free(tributaries);
    tributaries = NULL;
    return 0;
}

int main(void)
{
    int failed = 0;

    failed += test_openbound_outflow_removal_rate();
    failed += test_fischer_tributary_boundary_cell();

    if (failed > 0) {
        fprintf(stderr, "FAIL: %d test(s) failed\n", failed);
        return 1;
    }

    return 0;
}
