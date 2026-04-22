/**
 * @file test_calibration_target_summary.c
 * @brief Unit tests for Phase 1 calibration target summary helper.
 */

#include <stdio.h>
#include <string.h>

#include "variables.h"
#include "calibration_helpers.h"

static int fail(const char *msg)
{
    fprintf(stderr, "FAIL: %s\n", msg);
    return 1;
}

int main(void)
{
    CGEM_WQTargetSummaryInput targets[] = {
        // Stage 3: mix of chems + Phy total
        { 0,      3, CGEM_WQ_TARGET_KIND_CHEM,      O2  },
        { 3600,   3, CGEM_WQ_TARGET_KIND_CHEM,      NH4 },
        { 7200,   3, CGEM_WQ_TARGET_KIND_CHEM,      NH4 },
        { 86400,  3, CGEM_WQ_TARGET_KIND_CHEM,      NO3 },
        { 200000, 3, CGEM_WQ_TARGET_KIND_PHY_TOTAL, Phy1 },

        // Stage 2: simple two-point set (even count to exercise median convention)
        { 100,    2, CGEM_WQ_TARGET_KIND_CHEM,      SPM },
        { 200,    2, CGEM_WQ_TARGET_KIND_CHEM,      PO4 },
    };

    CGEM_WQTargetSummary sum;
    memset(&sum, 0, sizeof(sum));

    if (cgem_compute_wq_target_summary(targets, sizeof(targets)/sizeof(targets[0]), &sum) != 0) {
        return fail("cgem_compute_wq_target_summary returned error");
    }

    if (sum.total_targets != 7) return fail("total_targets should be 7");

    // Stage 3 checks
    if (sum.stage[3].total != 5) return fail("stage 3 total should be 5");
    if (sum.stage[3].by_chem[NH4] != 2) return fail("stage 3 NH4 count should be 2");
    if (sum.stage[3].by_chem[NO3] != 1) return fail("stage 3 NO3 count should be 1");
    if (sum.stage[3].by_chem[O2]  != 1) return fail("stage 3 O2 count should be 1");
    if (sum.stage[3].phy_total    != 1) return fail("stage 3 phy_total count should be 1");

    if (sum.stage[3].min_step_seconds != 0) return fail("stage 3 min_step_seconds should be 0");
    if (sum.stage[3].median_step_seconds != 7200) return fail("stage 3 median_step_seconds should be 7200");
    if (sum.stage[3].max_step_seconds != 200000) return fail("stage 3 max_step_seconds should be 200000");
    if (sum.stage[3].last_step_seconds != 200000) return fail("stage 3 last_step_seconds should equal max");

    // Stage 2 checks (median should be avg of {100,200} = 150)
    if (sum.stage[2].total != 2) return fail("stage 2 total should be 2");
    if (sum.stage[2].min_step_seconds != 100) return fail("stage 2 min_step_seconds should be 100");
    if (sum.stage[2].median_step_seconds != 150) return fail("stage 2 median_step_seconds should be 150");
    if (sum.stage[2].max_step_seconds != 200) return fail("stage 2 max_step_seconds should be 200");
    if (sum.stage[2].last_step_seconds != 200) return fail("stage 2 last_step_seconds should equal max");

    printf("PASS: calibration target summary helper\n");
    return 0;
}
