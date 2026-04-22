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

    M = 11;
    M1 = 10;
    M2 = 9;
    M3 = 11;

    LC1 = 25000.0;
    LC2 = 25000.0;
    index_2 = -1;
    segment_transition_width = 1;

    C_VDB = 0.05;
    D0_CORRECTION = 1.0;
    dispersion_abs_min = 1.0;
    dispersion_abs_max = 1e6;

    for (int i = 1; i <= M; ++i) {
        totalArea[i] = 1000.0;
        width[i] = 100.0;
        waterDepth[i] = 5.0;
        riverbed_depth[i] = 5.0;
        velocity[i] = 0.0;
        disp[i] = 0.0;
    }

    enable_upstream_discharge_smoothing = 0;
    upstream_discharge_scale = 1.0;

    // Minimal discharge forcing for Discharge_ups()
    static double discharge_time[2] = {0.0, 1.0};
    static double discharge_val[2] = {100.0, 100.0};
    forcingData[FORCING_DISCHARGE].time = discharge_time;
    forcingData[FORCING_DISCHARGE].data = discharge_val;
    forcingData[FORCING_DISCHARGE].dataSize = 2;
}

static int nearly_equal(double a, double b, double tol)
{
    return fabs(a - b) <= tol;
}

static void run_disp_profile(int mode, double *out)
{
    tributary_mass_injection_mode = (TributaryMassInjectionMode)mode;
    Dispcoef(0);
    for (int i = 1; i <= M; ++i) {
        out[i] = disp[i];
    }
}

int main(void)
{
    init_minimal_grid();

    tributaryEnabled = 1;
    numTributaries = 1;
    tributaries = (Tributary *)calloc(1, sizeof(Tributary));
    if (!tributaries) {
        fprintf(stderr, "FAIL: could not allocate tributaries\n");
        return 1;
    }

    tributaries[0].cellIndex = 7;
    static double trib_q_time[2] = {0.0, 1.0};
    static double trib_q_val[2] = {40.0, 40.0};
    tributaries[0].dischargeTime = trib_q_time;
    tributaries[0].discharge = trib_q_val;
    tributaries[0].dischargeDataSize = 2;

    // 1) Mode-independence parity: Dispcoef should use consistent freshwater magnitude
    // regardless of tributary mass injection mode.
    double disp_fischer[MAXM + 1] = {0.0};
    double disp_rutherford[MAXM + 1] = {0.0};

    run_disp_profile(TRIB_INJECTION_FISCHER, disp_fischer);
    run_disp_profile(TRIB_INJECTION_RUTHERFORD, disp_rutherford);

    for (int i = 2; i <= M; ++i) {
        if (!nearly_equal(disp_fischer[i], disp_rutherford[i], 1e-12)) {
            fprintf(stderr,
                    "FAIL: Dispcoef mode parity mismatch at i=%d (Fischer=%.15g Rutherford=%.15g)\n",
                    i, disp_fischer[i], disp_rutherford[i]);
            free(tributaries);
            tributaries = NULL;
            return 1;
        }
    }

    // 2) Freshwater sensitivity: larger tributary inflow should strengthen upstream
    // dispersion decay (smaller upstream D), not collapse to a mode-dependent floor.
    double disp_low_q[MAXM + 1] = {0.0};
    double disp_high_q[MAXM + 1] = {0.0};

    trib_q_val[0] = trib_q_val[1] = 0.0;
    run_disp_profile(TRIB_INJECTION_FISCHER, disp_low_q);

    trib_q_val[0] = trib_q_val[1] = 80.0;
    run_disp_profile(TRIB_INJECTION_FISCHER, disp_high_q);

    const double ratio_low_q = disp_low_q[M] / fmax(disp_low_q[1], 1e-12);
    const double ratio_high_q = disp_high_q[M] / fmax(disp_high_q[1], 1e-12);

    if (!(ratio_high_q < ratio_low_q)) {
        fprintf(stderr,
                "FAIL: expected stronger freshwater inflow to increase relative upstream decay (ratio_highQ=%.15g ratio_lowQ=%.15g)\\n",
                ratio_high_q, ratio_low_q);
        free(tributaries);
        tributaries = NULL;
        return 1;
    }

    free(tributaries);
    tributaries = NULL;
    return 0;
}
