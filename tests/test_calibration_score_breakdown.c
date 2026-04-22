/**
 * @file test_calibration_score_breakdown.c
 * @brief Unit tests for per-variable weighted scaled RMSE breakdown helper.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "variables.h"
#include "calibration_helpers.h"

static int fail(const char *msg)
{
    fprintf(stderr, "FAIL: %s\n", msg);
    return 1;
}

static int assert_near(const char *label, double got, double expect, double tol)
{
    if (!isfinite(got) || fabs(got - expect) > tol) {
        fprintf(stderr, "FAIL: %s got %.17g expect %.17g (tol %.3g)\n", label, got, expect, tol);
        return 1;
    }
    return 0;
}

int main(void)
{
    // Synthetic samples (already filtered for recorded+finite+stage match).
    // O2: d={1,2}, w={1,1}  -> sum_w=2, sum_w_sq=5, rmse=sqrt(5/2)
    // NH4: d={3}, w={2}     -> sum_w=2, sum_w_sq=18, rmse=3
    // Phy(total): d={4}, w={0.5} -> sum_w=0.5, sum_w_sq=8, rmse=4
    CGEM_ScoreSample s[] = {
        { CGEM_SCORE_SAMPLE_CHEM,      O2,  1.0, 1.0 },
        { CGEM_SCORE_SAMPLE_CHEM,      O2,  2.0, 1.0 },
        { CGEM_SCORE_SAMPLE_CHEM,      NH4, 3.0, 2.0 },
        { CGEM_SCORE_SAMPLE_PHY_TOTAL, Phy1, 4.0, 0.5 },
    };

    CGEM_ScoreBreakdown out;
    memset(&out, 0, sizeof(out));

    if (cgem_compute_score_breakdown(s, sizeof(s)/sizeof(s[0]), &out) != 0) {
        return fail("cgem_compute_score_breakdown returned error");
    }

    // Per-variable RMSE
    {
        const double rmse_o2 = cgem_score_agg_rmse_scaled(&out.by_chem[O2]);
        const double rmse_nh4 = cgem_score_agg_rmse_scaled(&out.by_chem[NH4]);
        const double rmse_phy = cgem_score_agg_rmse_scaled(&out.phy_total);

        if (assert_near("rmse(O2)", rmse_o2, sqrt(5.0/2.0), 1e-12)) return 1;
        if (assert_near("rmse(NH4)", rmse_nh4, 3.0, 1e-12)) return 1;
        if (assert_near("rmse(Phy(total))", rmse_phy, 4.0, 1e-12)) return 1;
    }

    // Contribution percentages should sum to ~100% (based on sum_w_sq).
    {
        const double total_sq = out.total_sum_w_sq;
        const double c_o2  = cgem_score_contribution_pct(out.by_chem[O2].sum_w_sq, total_sq);
        const double c_nh4 = cgem_score_contribution_pct(out.by_chem[NH4].sum_w_sq, total_sq);
        const double c_phy = cgem_score_contribution_pct(out.phy_total.sum_w_sq, total_sq);

        const double sum = c_o2 + c_nh4 + c_phy;
        if (assert_near("contrib_sum_pct", sum, 100.0, 1e-9)) return 1;
    }

    // Phy(total) bucket should be isolated from chem.
    if (out.phy_total.n != 1) return fail("phy_total.n should be 1");
    if (out.by_chem[Phy1].n != 0) return fail("chem bucket for Phy1 should remain 0 when using Phy(total) kind");

    printf("PASS: calibration score breakdown helper\n");
    return 0;
}
