/**
 * @file test_calibration_timebase_alignment.c
 * @brief Unit test: seasonal calibration targets use a post-warmup time base.
 */

#include <stdio.h>

#include "calibration_helpers.h"

static int fail(const char *msg)
{
    fprintf(stderr, "FAIL: %s\n", msg);
    return 1;
}

int main(void)
{
    const long day = 86400;
    const long warmup = 100 * day;
    const long target_step_seconds = 1 * day; // day 1 since post-warmup start

    // Must not capture during warmup, even if t_sim matches the target step seconds.
    if (cgem_should_capture_seasonal_target_at_sim_time(target_step_seconds, warmup, target_step_seconds)) {
        return fail("should not capture during warmup (t < WARMUP)");
    }

    // At the exact warmup boundary, t_rel=0, so should not capture a target at +1 day.
    if (cgem_should_capture_seasonal_target_at_sim_time(warmup, warmup, target_step_seconds)) {
        return fail("should not capture at t=WARMUP for target at +1 day");
    }

    // Must capture at simulation time t=WARMUP + target_step_seconds.
    if (!cgem_should_capture_seasonal_target_at_sim_time(warmup + target_step_seconds, warmup, target_step_seconds)) {
        return fail("should capture at t=WARMUP+target_step_seconds");
    }

    printf("PASS: calibration seasonal target timebase alignment\n");
    return 0;
}
