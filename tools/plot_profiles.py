#!/usr/bin/env python3
"""
plot_profiles.py — Comprehensive longitudinal water-quality profile validation.

Produces ONE multi-panel figure with all WQ parameters arranged in a grid.
Each panel shows simulated dry/wet seasonal mean ± std bands with optional
observed data overlay (scatter by season) and validation metrics (RMSE, R²).

Usage:
    python tools/plot_profiles.py                       # Model-only
    python tools/plot_profiles.py --with-obs             # Overlay observations
    python tools/plot_profiles.py --obs-file my_obs.csv  # Custom obs file
    python tools/plot_profiles.py --vars Sal O2 NH4      # Specific variables

Observation input can be either:
    1. Compact model-unit tables such as `obs_longitudinal_profiles.csv`
    2. Raw wide-format water-quality tables with common field/lab columns

Both are normalized internally before plotting.
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cgem_read
import validation_obs

DEFAULT_VARS = ["Sal", "O2", "NH4", "NO3", "PO4", "TOC", "SPM", "Phy1",
                "Si", "DIC", "AT", "pCO2", "PH"]
MIN_VALIDATION_DAYS_PER_SEASON = 14


def _stats(obs_km, obs_val, sim_km, sim_mean):
    """Compute phase3-style RMSE, Diff%, and R² by interpolating to obs locations."""
    obs_km = np.asarray(obs_km, dtype=np.float64)
    obs_val = np.asarray(obs_val, dtype=np.float64)
    sim_km = np.asarray(sim_km, dtype=np.float64)
    sim_mean = np.asarray(sim_mean, dtype=np.float64)
    mask = np.isfinite(obs_km) & np.isfinite(obs_val)
    if mask.sum() < 2:
        return np.nan, np.nan, np.nan
    obs_k, obs_v = obs_km[mask], obs_val[mask]
    sim_at_obs = np.interp(obs_k, sim_km, sim_mean)
    rmse = np.sqrt(np.mean((obs_v - sim_at_obs) ** 2))
    mean_obs = np.mean(obs_v)
    diff_pct = (100.0 * rmse / mean_obs) if abs(mean_obs) > 1e-12 else np.nan
    obs_centered = obs_v - obs_v.mean()
    sim_centered = sim_at_obs - sim_at_obs.mean()
    denom = np.sqrt(np.sum(obs_centered ** 2) * np.sum(sim_centered ** 2))
    r2 = ((obs_centered @ sim_centered) / denom) ** 2 if denom > 0 else np.nan
    return rmse, diff_pct, r2


def _validation_coverage_note(obs_df, sim_index):
    """Return a note when seasonal validation metrics would be misleading."""
    if obs_df is None or "Season" not in obs_df.columns:
        return None

    sim_dates = pd.DatetimeIndex(sim_index)
    if sim_dates.empty:
        return "Validation metrics skipped:\nno simulation dates found"

    sim_seasons = pd.Series(validation_obs.classify_validation_season(sim_dates), index=sim_dates)
    day_counts = {
        season_name: int(sim_dates[sim_seasons == season_name].normalize().nunique())
        for season_name in ("dry", "wet")
    }
    insufficient = [
        season_name for season_name, n_days in day_counts.items()
        if n_days < MIN_VALIDATION_DAYS_PER_SEASON
    ]
    if not insufficient:
        return None

    coverage = ", ".join(f"{sn}={day_counts[sn]} d" for sn in ("dry", "wet"))
    start = sim_dates.min().date().isoformat()
    end = sim_dates.max().date().isoformat()
    return (
        "Validation metrics skipped:\n"
        f"seasonal coverage too short\n"
        f"({start} → {end}; {coverage})"
    )


def plot(out_dir, obs_file=None, variables=None, save_dir=None):
    """Create comprehensive longitudinal profile figure."""
    out_dir = Path(out_dir)
    if save_dir is None:
        save_dir = out_dir
    fig_dir = Path(save_dir) / "figures"
    tbl_dir = Path(save_dir) / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir.mkdir(parents=True, exist_ok=True)

    available = cgem_read.list_variables(out_dir)
    if variables:
        plot_vars = [v for v in variables if v in available]
    else:
        plot_vars = [v for v in DEFAULT_VARS if v in available]
    if not plot_vars:
        print("  No model output variables found")
        return

    obs_df = None
    if obs_file is not None:
        obs_df, _ = validation_obs.load_standardized_observations(obs_file, kind="profiles")
    has_obs = obs_df is not None

    # --- Figure layout (phase3-inspired: 6.2 × 4.4 per panel) ---
    n = len(plot_vars)
    n_cols = min(3, n)
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(6.2 * n_cols, 4.4 * n_rows),
                             squeeze=False)

    all_metrics = []
    validation_notes = []

    for idx, var in enumerate(plot_vars):
        ax = axes[idx // n_cols][idx % n_cols]
        long_name, unit = cgem_read.get_variable_info(var)
        unit = validation_obs.get_unit_override(var) or unit
        ylabel = f"{long_name} ({unit})" if unit else long_name

        try:
            sim_df = cgem_read.load_output(var, out_dir)
        except FileNotFoundError:
            ax.text(0.5, 0.5, f"No {var} data file", ha="center", va="center",
                    transform=ax.transAxes, fontsize=11, color="gray")
            ax.set_title(ylabel)
            continue

        sim_scale = validation_obs.get_sim_scale(var)
        if sim_scale != 1.0:
            sim_df = sim_df.astype(float) * sim_scale

        model_km = sim_df.columns.astype(float).to_numpy(dtype=float)
        seasons = validation_obs.classify_validation_season(sim_df.index)
        validation_note = _validation_coverage_note(obs_df, sim_df.index)
        if validation_note:
            validation_notes.append(validation_note)

        # --- Simulated: overall mean ± std band ---
        overall_mean = sim_df.mean(axis=0).to_numpy(dtype=float)
        overall_std = sim_df.std(axis=0).to_numpy(dtype=float)
        ax.fill_between(model_km,
                        overall_mean - overall_std,
                        overall_mean + overall_std,
                        color="cornflowerblue", alpha=0.22)
        ax.plot(model_km, overall_mean, color="navy", lw=1.8, label="Simulated")

        # --- Simulated: dry/wet seasonal overlays (dashed) ---
        season_styles = {"dry": ("#d95f02", "--", "Sim dry"),
                         "wet": ("#1b9e77", "--", "Sim wet")}
        for sn, (color, ls, lbl) in season_styles.items():
            mask = seasons == sn
            if not mask.any():
                continue
            profile = sim_df.loc[mask].mean(axis=0).to_numpy(dtype=float)
            ax.plot(model_km, profile, color=color, lw=1.2, ls=ls,
                    alpha=0.85, label=lbl)

        # --- Observed data overlay ---
        if has_obs and var in obs_df.columns and "Distance_km" in obs_df.columns:
            obs_var = obs_df[["Distance_km", var]].copy()
            if "Season" in obs_df.columns:
                obs_var["Season"] = validation_obs.normalize_season_labels(obs_df["Season"])
            obs_var[var] = pd.to_numeric(obs_var[var], errors="coerce")
            obs_var = obs_var.dropna(subset=["Distance_km", var])

            # Scatter by season
            obs_styles = {"dry": ("#d95f02", "o", "Obs dry"),
                          "wet": ("#1b9e77", "s", "Obs wet")}
            for sn, (color, marker, lbl) in obs_styles.items():
                if "Season" in obs_var.columns:
                    sub = obs_var[obs_var["Season"] == sn]
                else:
                    sub = obs_var
                if sub.empty:
                    continue
                ax.scatter(sub["Distance_km"], sub[var], s=22, c=color,
                           marker=marker, label=lbl, zorder=5)

            # Unseasoned obs (if no Season column)
            if "Season" not in obs_var.columns and not obs_var.empty:
                ax.scatter(obs_var["Distance_km"], obs_var[var], s=18,
                           c="crimson", marker="x", label="Obs", zorder=5)

            # --- Metrics text box ---
            if validation_note:
                ax.text(0.02, 0.98, validation_note,
                        transform=ax.transAxes, ha="left", va="top",
                        fontsize=8.2,
                        bbox=dict(facecolor="white", edgecolor="0.7", alpha=0.9))
            else:
                rmse, diff_pct, r2 = _stats(
                    obs_var["Distance_km"].values,
                    obs_var[var].values,
                    model_km,
                    overall_mean,
                )
                if not np.isnan(rmse):
                    ax.text(0.02, 0.98,
                            f"RMSE: {rmse:.2f}\nDiff%: {diff_pct:.1f}\nR²: {r2:.3f}",
                            transform=ax.transAxes, ha="left", va="top",
                            fontsize=8.5,
                            bbox=dict(facecolor="white", edgecolor="0.7", alpha=0.85))
                    all_metrics.append({"Variable": var, "RMSE": rmse,
                                        "Diff_pct": diff_pct, "R2": r2})

        # --- Axes formatting ---
        ax.set_title(ylabel)
        ax.set_xlabel("Distance from mouth (km)")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(loc="best", fontsize=8.5, frameon=True, ncol=1)

    # Hide unused panels
    for idx in range(n, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)

    # --- Save figure ---
    try:
        plt.tight_layout(pad=2.0)
    except Exception:
        plt.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.08,
                            hspace=0.4, wspace=0.3)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fig_path = fig_dir / "fig_profiles.png"
    fig.savefig(fig_path, dpi=250, bbox_inches="tight")
    fig.savefig(fig_dir / f"fig_profiles_{ts}.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # --- Save metrics CSV ---
    metrics_df = pd.DataFrame(all_metrics)
    if not metrics_df.empty:
        metrics_df.to_csv(tbl_dir / "validation_profiles.csv",
                          index=False, float_format="%.4f")

    # --- Console summary ---
    print(f"\n{'='*65}")
    print("  LONGITUDINAL PROFILES")
    print(f"{'='*65}")
    if metrics_df.empty:
        if validation_notes:
            print("  Validation metrics skipped (insufficient seasonal coverage)")
            for note in sorted(set(validation_notes)):
                print(f"    - {note.replace(chr(10), '; ')}")
        else:
            status = "Obs loaded but no matching variables" if has_obs else "Model-only"
            print(f"  {status}")
    else:
        for _, row in metrics_df.iterrows():
            print(f"  {row['Variable']:12s}  RMSE={row['RMSE']:8.3f}  "
                  f"Diff={row['Diff_pct']:+7.1f}%  R²={row['R2']:.3f}")
    print(f"  Figure → {fig_path}")
    print(f"{'='*65}")


if __name__ == "__main__":
    DEFAULT_OBS = validation_obs.default_obs_path("profiles")
    parser = argparse.ArgumentParser(description="Longitudinal WQ profile validation")
    parser.add_argument("--output-dir", default="OUT", help="Model output directory")
    parser.add_argument("--with-obs", action="store_true",
                        help=f"Enable obs overlay using default file ({DEFAULT_OBS})")
    parser.add_argument("--obs-file", default=None,
                        help="Observation CSV (Distance_km, Season, Var1, ...)")
    parser.add_argument("--vars", nargs="*", default=None,
                        help="Variables to plot (default: all available)")
    parser.add_argument("--save-dir", default=None,
                        help="Output directory for figures/tables")
    args = parser.parse_args()
    obs = args.obs_file or (str(DEFAULT_OBS) if args.with_obs else None)
    plot(Path(args.output_dir), obs, args.vars,
         args.save_dir or args.output_dir)
