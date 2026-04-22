#ifndef CGEM_CALIBRATION_HELPERS_H
#define CGEM_CALIBRATION_HELPERS_H

/**
 * @file calibration_helpers.h
 * @brief Pure helper functions for calibration scoring, target summary, and
 *        timebase alignment.  Independent of NLopt and calibration iteration
 *        state — always compiled regardless of CGEM_ENABLE_CALIBRATION.
 */

#include <stdbool.h>
#include <stddef.h>
#include "variables.h"  /* Chem, CHEM_COUNT */

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================
 * Score breakdown (weighted scaled RMSE per variable)
 * ==================================================================== */

typedef enum {
    CGEM_SCORE_SAMPLE_CHEM = 0,
    CGEM_SCORE_SAMPLE_PHY_TOTAL = 1,
} CGEM_ScoreVarKind;

typedef struct {
    CGEM_ScoreVarKind kind;
    Chem chem;        /* valid when kind == CGEM_SCORE_SAMPLE_CHEM */
    double d;         /* scaled residual: (model - target) / scale */
    double weight;    /* target weight */
} CGEM_ScoreSample;

typedef struct {
    double sum_w;
    double sum_w_sq;
    size_t n;
} CGEM_ScoreAgg;

typedef struct {
    CGEM_ScoreAgg by_chem[CHEM_COUNT];
    CGEM_ScoreAgg phy_total;
    double total_sum_w;
    double total_sum_w_sq;
    size_t total_n;
} CGEM_ScoreBreakdown;

int    cgem_compute_score_breakdown(const CGEM_ScoreSample *samples, size_t n,
                                   CGEM_ScoreBreakdown *out);
double cgem_score_agg_rmse_scaled(const CGEM_ScoreAgg *a);
double cgem_score_contribution_pct(double sum_w_sq, double total_sum_w_sq);

/* ====================================================================
 * Target summary (per-stage counts + time statistics)
 * ==================================================================== */

typedef enum {
    CGEM_WQ_TARGET_KIND_CHEM = 0,
    CGEM_WQ_TARGET_KIND_PHY_TOTAL = 1,
} CGEM_WQTargetKind;

typedef struct {
    long target_step_seconds;
    int  stage;  /* 1..4 */
    CGEM_WQTargetKind kind;
    Chem chem;   /* valid if kind == CGEM_WQ_TARGET_KIND_CHEM */
} CGEM_WQTargetSummaryInput;

typedef struct {
    size_t total;
    size_t by_chem[CHEM_COUNT];
    size_t phy_total;
    long min_step_seconds;
    long median_step_seconds;
    long max_step_seconds;
    long last_step_seconds;
} CGEM_WQStageSummary;

typedef struct {
    size_t total_targets;
    CGEM_WQStageSummary stage[5];  /* index by stage 1..4 */
} CGEM_WQTargetSummary;

int cgem_compute_wq_target_summary(const CGEM_WQTargetSummaryInput *targets,
                                   size_t n, CGEM_WQTargetSummary *out);

/* ====================================================================
 * Timebase alignment (post-warmup capture predicate)
 * ==================================================================== */

bool cgem_should_capture_seasonal_target_at_sim_time(long t_sim_seconds,
                                                     long warmup_seconds,
                                                     long target_step_seconds);

#ifdef __cplusplus
}
#endif

#endif /* CGEM_CALIBRATION_HELPERS_H */
