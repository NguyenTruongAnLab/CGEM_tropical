#!/usr/bin/env python3
"""
plot_seasonal.py — Comprehensive seasonal time-series validation.

Produces ONE multi-panel figure: params × stations grid.
Each cell shows the tidally-filtered time series with min/max envelope,
optional observed data overlay, and validation metrics (RMSE, R²).

If model output is hourly, a Godin (1972) tidal filter is applied
(A24·A24·A25 running mean — removes M2, S2, K1, O1, N2, M4 constituents).

Usage:
    python tools/plot_seasonal.py                          # Model-only
    python tools/plot_seasonal.py --with-obs                # Overlay obs
    python tools/plot_seasonal.py --obs-file my_obs.csv     # Custom obs file
    python tools/plot_seasonal.py --vars Sal O2 NH4         # Specific vars
    python tools/plot_seasonal.py --stations 20 60 100 150  # Custom stations

Observation input can be either:
    1. Compact model-unit tables such as `obs_timeseries.csv`
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
import matplotlib.dates as mdates

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cgem_read
import validation_obs

DEFAULT_VARS = ["Sal", "O2", "NH4", "NO3", "PO4", "TOC", "SPM", "Phy1", "pCO2"]
DEFAULT_STATIONS_KM = [20, 60, 100, 140]
MIN_VALIDATION_DAYS = 14


def _godin_filter(series):
    """Godin (1972) tidal filter: A24·A24·A25 running mean."""
    s1 = series.rolling(window=24, center=True, min_periods=18).mean()
    s2 = s1.rolling(window=24, center=True, min_periods=18).mean()
    return s2.rolling(window=25, center=True, min_periods=19).mean()


def _detect_stations(obs_df):
    """Auto-detect station km from obs data, sorted by distance."""
    if obs_df is None or "Distance_km" not in obs_df.columns:
        return None
    stations = obs_df.groupby("Distance_km").size().reset_index(name="count")
    km_values = sorted(stations["Distance_km"].unique())
    # Limit to ~8 stations max for readability
    if len(km_values) > 8:
        indices = np.linspace(0, len(km_values) - 1, 8, dtype=int)
        km_values = [km_values[i] for i in indices]
    return km_values


def _station_label(km, obs_df):
    """Build station label (name + km) from obs data if available."""
    if obs_df is not None and "Station" in obs_df.columns and "Distance_km" in obs_df.columns:
        match = obs_df[obs_df["Distance_km"] == km]
        if not match.empty:
            name = match["Station"].iloc[0]
            return f"{name} (km {km:.0f})"
    return f"km {km:.0f}"


def _validation_window_note(sim_index):
    """Return a note when the simulation window is too short for robust metrics."""
    sim_dates = pd.DatetimeIndex(sim_index)
    if sim_dates.empty:
        return "Metrics skipped:\nno simulation dates found"
    n_days = int(sim_dates.normalize().nunique())
    if n_days >= MIN_VALIDATION_DAYS:
        return None
    start = sim_dates.min().date().isoformat()
    end = sim_dates.max().date().isoformat()
    return f"Metrics skipped:\nshort run ({start} → {end}; {n_days} d)"


def plot(out_dir, obs_file=None, variables=None, stations_km=None,
         save_dir=None):
    """Create comprehensive seasonal time-series figure."""
    out_dir = Path(out_dir)
    if save_dir is None:
        save_dir = out_dir
    fig_dir = Path(save_dir) / "figures"
    tbl_dir = Path(save_dir) / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir.mkdir(parents=True, exist_ok=True)

    # Load obs
    obs_df = None
    if obs_file is not None:
        obs_df, _ = validation_obs.load_standardized_observations(obs_file, kind="seasonal")
    has_obs = obs_df is not None

    # Determine variables
    available = cgem_read.list_variables(out_dir)
    if variables:
        plot_vars = [v for v in variables if v in available]
    else:
        plot_vars = [v for v in DEFAULT_VARS if v in available]
    if not plot_vars:
        print("  No model output variables found")
        return

    # Determine station locations
    if stations_km is None:
        auto = _detect_stations(obs_df)
        stations_km = auto if auto else DEFAULT_STATIONS_KM

    n_params = len(plot_vars)
    n_stations = len(stations_km)

    # --- Figure layout (phase4-inspired: params × stations grid) ---
    fig, axes = plt.subplots(n_params, n_stations,
                             figsize=(5 * n_stations, 2.8 * n_params),
                             squeeze=False)
    plt.subplots_adjust(hspace=0.35, wspace=0.28)

    all_metrics = []
    validation_note = None

    for pi, var in enumerate(plot_vars):
        long_name, unit = cgem_read.get_variable_info(var)
        unit = validation_obs.get_unit_override(var) or unit
        label = f"{long_name} ({unit})" if unit else long_name

        try:
            sim_df = cgem_read.load_output(var, out_dir)
        except FileNotFoundError:
            for si in range(n_stations):
                ax = axes[pi, si]
                ax.text(0.5, 0.5, f"No {var} data", ha="center", va="center",
                        transform=ax.transAxes, color="gray")
            continue

        sim_scale = validation_obs.get_sim_scale(var)
        if sim_scale != 1.0:
            sim_df = sim_df.astype(float) * sim_scale

        model_km = sim_df.columns.astype(float).values
        dates_index = sim_df.index
        validation_note = _validation_window_note(dates_index)

        # Detect hourly vs daily
        if len(dates_index) >= 2:
            median_dt = pd.Series(dates_index).diff().dropna().median()
            is_hourly = median_dt <= pd.Timedelta(hours=2)
        else:
            is_hourly = False

        for si, km in enumerate(stations_km):
            ax = axes[pi, si]
            # Find closest cell
            closest_idx = np.argmin(np.abs(model_km - km))
            closest_km = model_km[closest_idx]
            ts = sim_df.iloc[:, closest_idx]

            if is_hourly:
                # Godin-filtered mean + raw daily min/max envelope
                filtered = _godin_filter(ts)
                daily_groups = ts.groupby(ts.index.date)
                daily_min = daily_groups.min()
                daily_max = daily_groups.max()
                # Resample filtered to daily (noon)
                daily_mean = filtered.resample("D").mean().dropna()
                # Align all to same date range
                common_dates = sorted(set(daily_min.index) &
                                      set(daily_max.index) &
                                      set(daily_mean.index))
                if common_dates:
                    plot_dates = pd.to_datetime(common_dates)
                    y_min = daily_min.reindex(common_dates).values
                    y_max = daily_max.reindex(common_dates).values
                    y_mean = daily_mean.reindex(common_dates).values
                else:
                    plot_dates = dates_index
                    y_mean = filtered.values
                    y_min = y_max = y_mean
            else:
                # Daily data — use smoothed mean
                plot_dates = dates_index
                y_mean = ts.rolling(15, center=True, min_periods=8).mean().values
                y_min = ts.values
                y_max = ts.values

            # Plot envelope + mean
            ax.fill_between(plot_dates, y_min, y_max,
                            alpha=0.3, color="blue", linewidth=0)
            ax.plot(plot_dates, y_mean, color="blue", linewidth=1.5)

            # Obs overlay
            if has_obs and var in obs_df.columns and "Distance_km" in obs_df.columns:
                # Find obs near this station (within ±5 km)
                obs_near = obs_df[obs_df["Distance_km"] == km].dropna(subset=[var])
                if obs_near.empty:
                    obs_near = obs_df[
                        (obs_df["Distance_km"] - km).abs() <= 5.0
                    ].dropna(subset=[var])

                if not obs_near.empty and "Date" in obs_near.columns:
                    obs_vals = obs_near[var].astype(float).values
                    obs_dates = obs_near["Date"]

                    # Color by source if available
                    if "Source" in obs_near.columns:
                        for source, color, marker in [("CEM", "red", "o"),
                                                       ("CARE", "green", "o")]:
                            sub = obs_near[obs_near["Source"] == source]
                            if not sub.empty:
                                ax.scatter(sub["Date"], sub[var].astype(float),
                                           color=color, s=30, zorder=5,
                                           label=source, marker=marker)
                    else:
                        ax.scatter(obs_dates, obs_vals, color="red", s=25,
                                   zorder=5, marker="o", label="Obs")

                    # Compute RMSE, R² by matching obs to model at same dates
                    if validation_note:
                        ax.text(0.03, 0.95, validation_note,
                                transform=ax.transAxes, fontsize=8,
                                va="top",
                                bbox=dict(facecolor="white", alpha=0.8,
                                          edgecolor="none", pad=1))
                    else:
                        matched_sim = np.interp(
                            mdates.date2num(obs_dates),
                            mdates.date2num(plot_dates),
                            y_mean)
                        valid = np.isfinite(obs_vals) & np.isfinite(matched_sim)
                        if valid.sum() >= 2:
                            o, s = obs_vals[valid], matched_sim[valid]
                            rmse = np.sqrt(np.mean((o - s) ** 2))
                            ss_res = np.sum((o - s) ** 2)
                            ss_tot = np.sum((o - np.mean(o)) ** 2)
                            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
                            ax.text(0.03, 0.95,
                                    f"RMSE: {rmse:.2f}\nR²: {r2:.3f}",
                                    transform=ax.transAxes, fontsize=8,
                                    va="top",
                                    bbox=dict(facecolor="white", alpha=0.7,
                                              edgecolor="none", pad=1))
                            all_metrics.append({
                                "Variable": var, "Station_km": km,
                                "RMSE": rmse, "R2": r2, "N": int(valid.sum())})

            # Titles & labels
            if pi == 0:
                ax.set_title(_station_label(km, obs_df),
                             fontsize=10, fontweight="bold")
            if si == 0:
                ax.set_ylabel(label, fontsize=10, fontweight="bold")

            # X-axis formatting
            ax.set_xlim(plot_dates.min(), plot_dates.max())
            locator = mdates.AutoDateLocator(minticks=3, maxticks=8)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b-%Y"))
            if pi == n_params - 1:
                plt.setp(ax.xaxis.get_majorticklabels(),
                         rotation=45, fontsize=8, ha="right")
            else:
                ax.set_xticklabels([])

            ax.grid(True, linestyle=":", alpha=0.4)
            ax.set_ylim(bottom=0)

    # --- Save ---
    plt.tight_layout()
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fig_path = fig_dir / "fig_seasonal.png"
    fig.savefig(fig_path, dpi=250, bbox_inches="tight")
    fig.savefig(fig_dir / f"fig_seasonal_{ts}.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # Metrics CSV
    metrics_df = pd.DataFrame(all_metrics)
    if not metrics_df.empty:
        metrics_df.to_csv(tbl_dir / "validation_seasonal.csv",
                          index=False, float_format="%.4f")

    # Console summary
    print(f"\n{'='*65}")
    print("  SEASONAL TIME-SERIES VALIDATION")
    print(f"{'='*65}")
    if metrics_df.empty:
        if validation_note:
            print(f"  {validation_note.replace(chr(10), '; ')}")
        else:
            status = "Obs loaded but no matching data" if has_obs else "Model-only"
            print(f"  {status}")
    else:
        for _, row in metrics_df.iterrows():
            print(f"  {row['Variable']:12s}  km={row['Station_km']:5.0f}  "
                  f"N={row['N']:3.0f}  RMSE={row['RMSE']:8.3f}  "
                  f"R²={row['R2']:.3f}")
    print(f"  Figure → {fig_path}")
    print(f"{'='*65}")


if __name__ == "__main__":
    DEFAULT_OBS = validation_obs.default_obs_path("seasonal")
    parser = argparse.ArgumentParser(description="Seasonal time-series validation")
    parser.add_argument("--output-dir", default="OUT", help="Model output directory")
    parser.add_argument("--with-obs", action="store_true",
                        help=f"Enable obs overlay using default file ({DEFAULT_OBS})")
    parser.add_argument("--obs-file", default=None,
                        help="Observation CSV (Date, Station, Distance_km, ...)")
    parser.add_argument("--vars", nargs="*", default=None,
                        help="Variables to plot (default: Sal O2 NH4 NO3 PO4 TOC SPM Phy1 pCO2)")
    parser.add_argument("--stations", nargs="*", type=float, default=None,
                        help="Station distances from mouth in km")
    parser.add_argument("--save-dir", default=None,
                        help="Output directory for figures/tables")
    args = parser.parse_args()
    obs = args.obs_file or (str(DEFAULT_OBS) if args.with_obs else None)
    plot(Path(args.output_dir), obs, args.vars, args.stations,
         args.save_dir or args.output_dir)
