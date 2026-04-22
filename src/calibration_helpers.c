/**
 * @file calibration_helpers.c
 * @brief Pure helper functions for calibration scoring, target summary, and
 *        timebase alignment.  Always compiled (no NLopt dependency).
 */

#include "calibration_helpers.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ====================================================================
 * Score breakdown
 * ==================================================================== */

static void score_agg_add(CGEM_ScoreAgg *a, double d, double w)
{
    if (!a) return;
    a->sum_w   += w;
    a->sum_w_sq += w * d * d;
    a->n       += 1;
}

int cgem_compute_score_breakdown(const CGEM_ScoreSample *samples, size_t n,
                                 CGEM_ScoreBreakdown *out)
{
    if (!out) return -1;
    memset(out, 0, sizeof(*out));
    if (!samples && n > 0) return -1;

    for (size_t i = 0; i < n; ++i) {
        const CGEM_ScoreSample *s = &samples[i];
        if (!isfinite(s->d) || !isfinite(s->weight) || !(s->weight > 0.0))
            continue;

        out->total_sum_w    += s->weight;
        out->total_sum_w_sq += s->weight * s->d * s->d;
        out->total_n        += 1;

        if (s->kind == CGEM_SCORE_SAMPLE_PHY_TOTAL) {
            score_agg_add(&out->phy_total, s->d, s->weight);
        } else {
            const int c = (int)s->chem;
            if (c >= 0 && c < CHEM_COUNT)
                score_agg_add(&out->by_chem[c], s->d, s->weight);
        }
    }
    return 0;
}

double cgem_score_agg_rmse_scaled(const CGEM_ScoreAgg *a)
{
    if (!a) return NAN;
    if (!(a->sum_w > 0.0)) return NAN;
    if (!isfinite(a->sum_w_sq)) return NAN;
    return sqrt(a->sum_w_sq / a->sum_w);
}

double cgem_score_contribution_pct(double sum_w_sq, double total_sum_w_sq)
{
    if (!(total_sum_w_sq > 0.0) || !isfinite(total_sum_w_sq) ||
        !isfinite(sum_w_sq) || sum_w_sq < 0.0)
        return 0.0;
    return 100.0 * (sum_w_sq / total_sum_w_sq);
}

/* ====================================================================
 * Target summary
 * ==================================================================== */

static int compare_long_asc(const void *a, const void *b)
{
    const long va = *(const long *)a;
    const long vb = *(const long *)b;
    if (va < vb) return -1;
    if (va > vb) return  1;
    return 0;
}

int cgem_compute_wq_target_summary(const CGEM_WQTargetSummaryInput *targets,
                                   size_t n, CGEM_WQTargetSummary *out)
{
    if (!out) return -1;
    memset(out, 0, sizeof(*out));
    out->total_targets = n;
    if (!targets && n > 0) return -1;

    size_t stage_counts[5] = {0};
    for (size_t i = 0; i < n; ++i) {
        const int st = targets[i].stage;
        if (st >= 1 && st <= 4) stage_counts[st]++;
    }

    long *stage_steps[5] = {0};
    size_t stage_pos[5]  = {0};
    for (int st = 1; st <= 4; ++st) {
        if (stage_counts[st] > 0) {
            stage_steps[st] = (long *)malloc(stage_counts[st] * sizeof(long));
            if (!stage_steps[st]) {
                for (int s2 = 1; s2 <= 4; ++s2) free(stage_steps[s2]);
                return -1;
            }
        }
    }

    for (size_t i = 0; i < n; ++i) {
        const CGEM_WQTargetSummaryInput *t = &targets[i];
        const int st = t->stage;
        if (st < 1 || st > 4) continue;

        out->stage[st].total++;
        if (t->kind == CGEM_WQ_TARGET_KIND_PHY_TOTAL)
            out->stage[st].phy_total++;
        else {
            const int c = (int)t->chem;
            if (c >= 0 && c < CHEM_COUNT) out->stage[st].by_chem[c]++;
        }

        if (stage_steps[st] && stage_pos[st] < stage_counts[st])
            stage_steps[st][stage_pos[st]++] = t->target_step_seconds;
    }

    for (int st = 1; st <= 4; ++st) {
        CGEM_WQStageSummary *s = &out->stage[st];
        const size_t k = s->total;
        if (k == 0 || !stage_steps[st]) {
            s->min_step_seconds = s->median_step_seconds = 0;
            s->max_step_seconds = s->last_step_seconds   = 0;
            continue;
        }
        qsort(stage_steps[st], k, sizeof(long), compare_long_asc);
        s->min_step_seconds  = stage_steps[st][0];
        s->max_step_seconds  = stage_steps[st][k - 1];
        s->last_step_seconds = s->max_step_seconds;
        if (k % 2 == 1)
            s->median_step_seconds = stage_steps[st][k / 2];
        else {
            const long a = stage_steps[st][(k / 2) - 1];
            const long b = stage_steps[st][k / 2];
            s->median_step_seconds = (long)llround(0.5 * ((double)a + (double)b));
        }
    }

    for (int st = 1; st <= 4; ++st) free(stage_steps[st]);
    return 0;
}

/* ====================================================================
 * Timebase alignment
 * ==================================================================== */

bool cgem_should_capture_seasonal_target_at_sim_time(long t_sim_seconds,
                                                     long warmup_seconds,
                                                     long target_step_seconds)
{
    if (warmup_seconds < 0) warmup_seconds = 0;
    if (target_step_seconds < 0) return false;
    if (t_sim_seconds < warmup_seconds) return false;
    return (t_sim_seconds - warmup_seconds) >= target_step_seconds;
}
