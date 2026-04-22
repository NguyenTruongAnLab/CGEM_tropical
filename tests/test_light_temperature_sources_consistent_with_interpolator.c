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
    WARMUP = 86400;

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

    const int times[] = {
        (int)WARMUP - 10000,
        (int)WARMUP - 1,
        (int)WARMUP,
        (int)WARMUP + 1,
        (int)WARMUP + 10000,
        (int)WARMUP + 10 * 86400
    };

    for (size_t k = 0; k < sizeof(times) / sizeof(times[0]); ++k) {
        const int t = times[k];

        const double T_interp = interpolateInputData(t,
                                                    forcingData[FORCING_TEMPERATURE].time,
                                                    forcingData[FORCING_TEMPERATURE].data,
                                                    forcingData[FORCING_TEMPERATURE].dataSize,
                                                    false,
                                                    "Temperature",
                                                    20.0);
        const double T_func = waterT(t);
        if (!nearly_equal(T_func, T_interp, 1e-12)) {
            fprintf(stderr, "FAIL: waterT != interpolateInputData at t=%d (got %.15g expected %.15g)\n",
                    t, T_func, T_interp);
            return 1;
        }

        const double I_interp = interpolateInputData(t,
                                                    forcingData[FORCING_LIGHT].time,
                                                    forcingData[FORCING_LIGHT].data,
                                                    forcingData[FORCING_LIGHT].dataSize,
                                                    true,
                                                    "Light",
                                                    0.0);
        const double I_func = I0(t);
        if (!nearly_equal(I_func, I_interp, 1e-12)) {
            fprintf(stderr, "FAIL: I0 != interpolateInputData at t=%d (got %.15g expected %.15g)\n",
                    t, I_func, I_interp);
            return 1;
        }
    }

    return 0;
}
