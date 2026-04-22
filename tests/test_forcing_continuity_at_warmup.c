#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "define.h"
#include "variables.h"

// Implemented in src/bcforcing.c
extern double interpolateInputData(int t, double* timeArray, double* dataArray, int dataSize,
                                  bool isHourly, const char* dataName, double defaultValue);

// Implemented in src/biogeout.c
extern double waterT(int t);
extern double I0(int t);

static int nearly_equal(double a, double b, double tol)
{
    return fabs(a - b) <= tol;
}

int main(void)
{
    // 1-day warmup for a focused, deterministic unit test.
    WARMUP = 86400;

    // Synthetic daily temperature (days) and hourly light (hours) datasets.
    // Chosen to be non-equal at endpoints so any non-periodic wrap would cause a jump.
    static double temp_time_days[3] = {0.0, 1.0, 2.0};
    static double temp_vals[3]      = {20.0, 21.0, 22.0};

    static double light_time_hours[3] = {0.0, 1.0, 2.0};
    static double light_vals[3]       = {100.0, 200.0, 300.0};

    forcingData[FORCING_TEMPERATURE].time = temp_time_days;
    forcingData[FORCING_TEMPERATURE].data = temp_vals;
    forcingData[FORCING_TEMPERATURE].dataSize = 3;

    forcingData[FORCING_LIGHT].time = light_time_hours;
    forcingData[FORCING_LIGHT].data = light_vals;
    forcingData[FORCING_LIGHT].dataSize = 3;

    const int t_before = (int)WARMUP - 1;
    const int t_at = (int)WARMUP;
    const int t_after = (int)WARMUP + 1;

    // Temperature continuity at warmup boundary.
    const double T_before = waterT(t_before);
    const double T_at = waterT(t_at);
    const double T_after = waterT(t_after);

    // Light continuity at warmup boundary.
    const double I_before = I0(t_before);
    const double I_at = I0(t_at);
    const double I_after = I0(t_after);

    // With periodic (circular) interpolation, the discontinuity at the wrap point should be O(dt).
    // Use conservative tolerances (dt = 1s).
    if (!nearly_equal(T_before, T_at, 1e-2)) {
        fprintf(stderr, "FAIL: waterT discontinuity at warmup (before=%.12g at=%.12g)\n", T_before, T_at);
        return 1;
    }
    if (!nearly_equal(I_before, I_at, 1.0)) {
        fprintf(stderr, "FAIL: I0 discontinuity at warmup (before=%.12g at=%.12g)\n", I_before, I_at);
        return 1;
    }

    // Sanity: values should remain finite and non-negative (light).
    if (!isfinite(T_before) || !isfinite(T_at) || !isfinite(T_after)) {
        fprintf(stderr, "FAIL: non-finite waterT near warmup\n");
        return 1;
    }
    if (!isfinite(I_before) || !isfinite(I_at) || !isfinite(I_after)) {
        fprintf(stderr, "FAIL: non-finite I0 near warmup\n");
        return 1;
    }
    if (I_before < 0.0 || I_at < 0.0 || I_after < 0.0) {
        fprintf(stderr, "FAIL: negative I0 near warmup (before=%.12g at=%.12g after=%.12g)\n", I_before, I_at, I_after);
        return 1;
    }

    return 0;
}
