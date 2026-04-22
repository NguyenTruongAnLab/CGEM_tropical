#!/usr/bin/env python3
"""
plot_hydrodynamics.py — Comprehensive estuarine transport diagnostic figure.

Produces ONE multi-panel figure (3 rows × 5 cols) covering geometry,
hydrodynamics, and transport diagnostics following the phase2 convention:

  Row 1: Riverbed depth, width, Chezy/friction, dispersion, mass balance
  Row 2: Water depth time series, water level contour (last 5 days),
         seasonal tidal range, particle trajectories, residence time
  Row 3: Salinity time series, seasonal salinity profile,
         salinity contour (2 tidal cycles), velocity + displacement,
         Eulerian vs Lagrangian salinity

With --with-obs, overlays observed tidal range and salinity data.

Usage:
    python tools/plot_hydrodynamics.py
    python tools/plot_hydrodynamics.py --with-obs
    python tools/plot_hydrodynamics.py --obs-file path/to/obs.csv

Observation CSV format (all values in model units):
    Distance_km, Tidal_range_m[, Salinity, Season, Date]
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
import matplotlib.gridspec as gridspec

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cgem_read


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rmse(obs, sim):
    return np.sqrt(np.nanmean((obs - sim) ** 2))


def _pearson_r(obs, sim):
    mask = np.isfinite(obs) & np.isfinite(sim)
    if mask.sum() < 3:
        return np.nan
    return np.corrcoef(obs[mask], sim[mask])[0, 1]


def _willmott_skill(obs, sim):
    mean_obs = np.nanmean(obs)
    num = np.nansum((obs - sim) ** 2)
    den = np.nansum((np.abs(sim - mean_obs) + np.abs(obs - mean_obs)) ** 2)
    return 1.0 - num / den if den > 0 else np.nan


def _tidal_range_per_cycle(depth_arr, axis_time, tidal_period_steps):
    """Compute tidal range statistics per M2 cycle."""
    n_steps = depth_arr.shape[0]
    n_cycles = n_steps // tidal_period_steps
    if n_cycles < 2:
        single = np.nanpercentile(depth_arr, 95, axis=0) - np.nanpercentile(depth_arr, 5, axis=0)
        return single, single, single
    ranges = []
    for i in range(n_cycles):
        s = depth_arr[i * tidal_period_steps:(i + 1) * tidal_period_steps]
        ranges.append(np.nanpercentile(s, 95, axis=0) -
                      np.nanpercentile(s, 5, axis=0))
    arr = np.array(ranges)
    return np.nanmean(arr, axis=0), np.nanpercentile(arr, 5, axis=0), np.nanpercentile(arr, 95, axis=0)


def _read_dt_seconds(out_dir):
    """Read model timestep from params.txt."""
    params = cgem_read.load_params()
    for key in ("DT", "dt", "Dt"):
        if key in params:
            try:
                return float(params[key])
            except ValueError:
                pass
    # Fallback: infer from output timestamps
    return 300.0


def _get_season_masks(dates):
    """Return (dry_mask, wet_mask) boolean arrays."""
    seasons = cgem_read.classify_season(dates)
    return seasons == "dry", seasons == "wet"


def _text_or_empty(ax, msg):
    """Show placeholder text on empty panel."""
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=11, color="gray",
            transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])


# ---------------------------------------------------------------------------
# Main plotting
# ---------------------------------------------------------------------------

def plot(out_dir, obs_file=None, save_dir=None):
    """Create comprehensive hydrodynamic diagnostic figure."""
    out_dir = Path(out_dir)
    if save_dir is None:
        save_dir = out_dir
    fig_dir = Path(save_dir) / "figures"
    tbl_dir = Path(save_dir) / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir.mkdir(parents=True, exist_ok=True)

    # Load obs
    obs_df = None
    if obs_file and Path(obs_file).exists():
        obs_df = pd.read_csv(obs_file)

    # Load geometry
    geom = cgem_read.load_geometry(out_dir)
    dist_km = geom["Distance_km"].values
    length_km = dist_km[-1] if len(dist_km) > 0 else 200.0

    # Load model outputs
    available = cgem_read.list_variables(out_dir)
    depth_df = sal_df = vel_df = disp_df = None
    dt_sec = _read_dt_seconds(out_dir)
    tidal_period_steps = max(1, int(round(44712.0 / dt_sec)))

    for var, name in [("waterDepth", "depth_df"), ("Sal", "sal_df"),
                      ("velocity", "vel_df"), ("disp", "disp_df")]:
        if var in available:
            try:
                locals()[name]  # just to suppress lint
            except Exception:
                pass
    try:
        depth_df = cgem_read.load_output("waterDepth", out_dir)
    except FileNotFoundError:
        pass
    try:
        sal_df = cgem_read.load_output("Sal", out_dir)
    except FileNotFoundError:
        pass
    try:
        vel_df = cgem_read.load_output("velocity", out_dir)
    except FileNotFoundError:
        pass
    try:
        disp_df = cgem_read.load_output("disp", out_dir)
    except FileNotFoundError:
        pass

    # ----- Create figure: 3 rows × 5 cols -----
    fig = plt.figure(figsize=(25, 15))
    gs = gridspec.GridSpec(3, 5, hspace=0.5, wspace=0.4)

    # ===== ROW 1: Geometry + Dispersion + Mass Balance =====
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[0, 3])
    ax5 = fig.add_subplot(gs[0, 4])

    # R1C1: Riverbed depth
    depth_col = "RiverbedDepth" if "RiverbedDepth" in geom.columns else "Depth"
    if depth_col in geom.columns:
        bd = geom[depth_col].values
        ax1.fill_between(dist_km, bd, 0, alpha=0.3, color="saddlebrown", label="Riverbed")
        ax1.plot(dist_km, bd, "k-", lw=1.2)
        ax1.axhline(0, color="steelblue", lw=0.8, ls="--", label="Water surface")
        ax1.invert_yaxis()
    ax1.set(title="Riverbed Depth", xlabel="Dist (km)", ylabel="Elevation (m)")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # R1C2: Width
    if "Width" in geom.columns:
        w = geom["Width"].values
        ax2.fill_between(dist_km, 0, w / 1000, alpha=0.2, color="steelblue")
        ax2.plot(dist_km, w / 1000, "b-", lw=1.2)
    ax2.set(title="Channel Width", xlabel="Dist (km)", ylabel="Width (km)")
    ax2.grid(True, alpha=0.3)

    # R1C3: Chezy or cross-section
    if "Chezy" in geom.columns:
        ax3.plot(dist_km, geom["Chezy"].values, "g-", lw=1.2)
        ax3.set(title="Chezy Coefficient", xlabel="Dist (km)", ylabel="C (m^½/s)")
    elif "FRIC" in geom.columns:
        ax3.plot(dist_km, geom["FRIC"].values, "g-", lw=1.2)
        ax3.set(title="Friction", xlabel="Dist (km)")
    else:
        if "Width" in geom.columns and depth_col in geom.columns:
            area = np.abs(geom[depth_col].values) * geom["Width"].values
            ax3.fill_between(dist_km, 0, area, alpha=0.2, color="green")
            ax3.plot(dist_km, area, "g-", lw=1.2)
            ax3.set(title="Cross-Section Area", xlabel="Dist (km)", ylabel="Area (m²)")
    ax3.grid(True, alpha=0.3)

    # R1C4: Dispersion
    if disp_df is not None and not disp_df.empty:
        disp_mean = disp_df.mean(axis=0).values
        disp_p5 = disp_df.quantile(0.05, axis=0).values
        disp_p95 = disp_df.quantile(0.95, axis=0).values
        disp_km = disp_df.columns.astype(float).values
        ax4.plot(disp_km, disp_mean, "b-", lw=1.2, label="Mean")
        ax4.fill_between(disp_km, disp_p5, disp_p95, color="b", alpha=0.2, label="5-95%")
        ax4.set(title="Dispersion Profile", xlabel="Dist (km)", ylabel="D (m²/s)")
        ax4.legend(fontsize=8)
    else:
        _text_or_empty(ax4, "No Dispersion Data")
    ax4.grid(True, alpha=0.3)

    # R1C5: Placeholder for mass balance or domain info
    info_text = f"Domain: {length_km:.0f} km\nCells: {len(dist_km)}\nΔx: {dist_km[1]-dist_km[0]:.1f} km" if len(dist_km) > 1 else "No geometry"
    params = cgem_read.load_params()
    for key in ("WARMUP", "MAXT"):
        if key in params:
            info_text += f"\n{key}: {params[key]}"
    ax5.text(0.1, 0.9, info_text, transform=ax5.transAxes, fontsize=11,
             va="top", family="monospace",
             bbox=dict(facecolor="lightyellow", edgecolor="gray", alpha=0.8))
    ax5.set_title("Model Domain")
    ax5.set_xticks([])
    ax5.set_yticks([])

    # ===== ROW 2: Depth + Water Level + Tidal Range + Trajectories + Residence Time =====
    ax_r2 = [fig.add_subplot(gs[1, i]) for i in range(5)]

    if depth_df is not None and not depth_df.empty:
        depth_km = depth_df.columns.astype(float).values
        depth_arr = depth_df.values
        depth_dates = depth_df.index

        # R2C1: Water depth time series at selected cells
        n_cells = depth_arr.shape[1]
        indices = np.linspace(0, n_cells - 1, min(4, n_cells), dtype=int)
        time_days = (depth_dates - depth_dates[0]).total_seconds() / 86400.0
        for i in indices:
            ax_r2[0].plot(time_days, depth_arr[:, i], label=f"{depth_km[i]:.0f} km")
        ax_r2[0].set(title="Water Depth Time Series", xlabel="Time (days)", ylabel="Depth (m)")
        ax_r2[0].legend(fontsize=8)
        ax_r2[0].grid(True, alpha=0.3)

        # R2C2: Water level contour (last 5 days)
        last_5d_mask = time_days >= (time_days.values[-1] - 5.0) if hasattr(time_days, 'values') else time_days >= (time_days[-1] - 5.0)
        if depth_col in geom.columns:
            bed = geom[depth_col].values
            wl = depth_arr[last_5d_mask, :] - bed[np.newaxis, :len(depth_arr[0])]
            wl_max = np.nanpercentile(np.abs(wl), 98)
            if wl_max > 0:
                c = ax_r2[1].contourf(
                    *np.meshgrid(depth_km, time_days[last_5d_mask]),
                    wl, levels=np.linspace(-wl_max, wl_max, 21), cmap="RdBu_r")
                fig.colorbar(c, ax=ax_r2[1], label="Level (m)")
        ax_r2[1].set(title="Water Level (Last 5 Days)", xlabel="Dist (km)", ylabel="Day")
        ax_r2[1].grid(True, alpha=0.3)

        # R2C3: Seasonal tidal range
        dry_mask, wet_mask = _get_season_masks(depth_dates)
        r2_dry_txt, r2_wet_txt = "", ""
        for mask, color, sn_label in [(dry_mask, "r", "Dry"), (wet_mask, "b", "Wet")]:
            if mask.any():
                mean_r, p5_r, p95_r = _tidal_range_per_cycle(
                    depth_arr[mask], None, tidal_period_steps)
                ax_r2[2].plot(depth_km, mean_r, f"{color[0]}-", label=f"Sim {sn_label}")
                ax_r2[2].fill_between(depth_km, p5_r, p95_r, color=color, alpha=0.1)

        # Obs tidal range overlay
        if obs_df is not None and "Tidal_range_m" in obs_df.columns:
            ax_r2[2].scatter(obs_df["Distance_km"], obs_df["Tidal_range_m"],
                             c="k", marker="o", s=20, alpha=0.7, label="Obs", zorder=5)
            # Compute R² against overall mean
            overall_mean, _, _ = _tidal_range_per_cycle(depth_arr, None, tidal_period_steps)
            sim_at_obs = np.interp(obs_df["Distance_km"].values, depth_km, overall_mean)
            obs_vals = obs_df["Tidal_range_m"].values
            r = _pearson_r(obs_vals, sim_at_obs)
            if not np.isnan(r):
                ax_r2[2].set_title(f"Seasonal Tidal Range ($R^2$={r**2:.2f})")
            else:
                ax_r2[2].set_title("Seasonal Tidal Range")
        else:
            ax_r2[2].set_title("Seasonal Tidal Range")
        ax_r2[2].set(xlabel="Dist (km)", ylabel="Range (m)")
        ax_r2[2].legend(fontsize=8)
        ax_r2[2].grid(True, alpha=0.3)

        # R2C4: Particle trajectories
        if vel_df is not None and not vel_df.empty:
            vel_arr = vel_df.values
            vel_km = vel_df.columns.astype(float).values
            vel_time_days = (vel_df.index - vel_df.index[0]).total_seconds() / 86400.0
            dt_days = dt_sec / 86400.0
            start_locs = np.linspace(0, length_km, 5)
            positions = np.zeros((len(vel_time_days), len(start_locs)))
            positions[0, :] = start_locs
            for t in range(len(vel_time_days) - 1):
                vel_at_pos = np.array([np.interp(p, vel_km, vel_arr[t]) for p in positions[t]])
                positions[t + 1] = np.clip(
                    positions[t] + vel_at_pos * 86400.0 / 1000.0 * dt_days, 0, length_km)
            for i in range(len(start_locs)):
                ax_r2[3].plot(vel_time_days, positions[:, i],
                              label=f"Start: {start_locs[i]:.0f} km")
            ax_r2[3].set(title="Particle Trajectories", xlabel="Time (days)", ylabel="Position (km)")
            ax_r2[3].legend(fontsize=7)
            ax_r2[3].grid(True, ls=":", alpha=0.4)
        else:
            _text_or_empty(ax_r2[3], "No Velocity Data")

        # R2C5: Residence time estimates
        if vel_df is not None and not vel_df.empty:
            residual_vel = np.nanmean(vel_arr, axis=0)  # m/s
            locs = np.linspace(0, length_km, 5)
            rt = []
            for loc in locs:
                rv = np.interp(loc, vel_km, residual_vel)
                rt_days = ((length_km - loc) * 1000.0) / (abs(rv) * 86400.0 + 1e-9)
                rt.append(rt_days)
            ax_r2[4].bar([f"{l:.0f}" for l in locs], rt, color="gray", edgecolor="k")
            ax_r2[4].set(title="Est. Residence Time", xlabel="Start km", ylabel="Time (days)")
            ax_r2[4].grid(True, axis="y", alpha=0.3)
        else:
            _text_or_empty(ax_r2[4], "No Velocity Data")
    else:
        for ax in ax_r2:
            _text_or_empty(ax, "No Water Depth Data")

    # ===== ROW 3: Salinity diagnostics =====
    ax_r3 = [fig.add_subplot(gs[2, i]) for i in range(5)]

    if sal_df is not None and not sal_df.empty:
        sal_arr = sal_df.values
        sal_km = sal_df.columns.astype(float).values
        sal_dates = sal_df.index
        sal_time_days = (sal_dates - sal_dates[0]).total_seconds() / 86400.0

        # R3C1: Salinity time series
        mean_profile = np.nanmean(sal_arr, axis=0)
        max_sal = np.nanmax(mean_profile) if mean_profile.size else 0
        if max_sal > 0:
            sal_threshold = max_sal * 0.05
            candidates = np.where(mean_profile >= sal_threshold)[0]
        else:
            candidates = np.arange(min(6, sal_arr.shape[1]))
        n_series = min(6, len(candidates))
        if n_series > 0:
            sel = np.linspace(0, len(candidates) - 1, n_series, dtype=int)
            for s in sel:
                idx = candidates[s]
                ax_r3[0].plot(sal_time_days, sal_arr[:, idx],
                              label=f"{sal_km[idx]:.0f} km")
        ax_r3[0].set(title="Salinity Time Series", xlabel="Time (days)", ylabel="Salinity (PSU)")
        ax_r3[0].legend(fontsize=7)
        ax_r3[0].grid(True, alpha=0.3)

        # R3C2: Seasonal salinity profile
        dry_mask, wet_mask = _get_season_masks(sal_dates)
        r2_dry, r2_wet = "", ""
        for mask, color, sn in [(dry_mask, "r", "Dry"), (wet_mask, "b", "Wet")]:
            if mask.any():
                s_mean = np.nanmean(sal_arr[mask], axis=0)
                s_p5 = np.nanpercentile(sal_arr[mask], 5, axis=0)
                s_p95 = np.nanpercentile(sal_arr[mask], 95, axis=0)
                ax_r3[1].plot(sal_km, s_mean, f"{color[0]}-", label=f"Sim {sn}")
                ax_r3[1].fill_between(sal_km, s_p5, s_p95, color=color, alpha=0.1)

        # Obs salinity overlay
        if obs_df is not None and "Salinity" in obs_df.columns:
            if "Season" in obs_df.columns:
                for sn, marker, color in [("dry", "x", "r"), ("wet", "+", "b")]:
                    sub = obs_df[obs_df["Season"] == sn].dropna(subset=["Salinity"])
                    if not sub.empty:
                        ax_r3[1].scatter(sub["Distance_km"], sub["Salinity"],
                                         c=color, marker=marker, s=30, alpha=0.9,
                                         label=f"Obs {sn.capitalize()}")
            else:
                obs_s = obs_df.dropna(subset=["Salinity"])
                if not obs_s.empty:
                    ax_r3[1].scatter(obs_s["Distance_km"], obs_s["Salinity"],
                                     c="k", marker="x", s=30, label="Obs")
        ax_r3[1].set(title="Seasonal Salinity Profile", xlabel="Dist (km)")
        ax_r3[1].legend(fontsize=8)
        ax_r3[1].grid(True, alpha=0.3)

        # R3C3: Salinity contour (last 2 tidal cycles)
        n_2cycles = min(2 * tidal_period_steps, sal_arr.shape[0])
        sal_slice = sal_arr[-n_2cycles:, :]
        time_hours = np.arange(n_2cycles) * dt_sec / 3600.0
        try:
            c = ax_r3[2].contourf(
                *np.meshgrid(sal_km, time_hours), sal_slice,
                cmap="magma", levels=20)
            fig.colorbar(c, ax=ax_r3[2], label="Salinity (PSU)")
        except Exception:
            pass
        ax_r3[2].set(title="Salinity (2 Tidal Cycles)", xlabel="Dist (km)", ylabel="Time (h)")

        # R3C4: Velocity + displacement at midpoint
        if vel_df is not None and not vel_df.empty:
            vel_arr2 = vel_df.values
            vel_km2 = vel_df.columns.astype(float).values
            mid_km = length_km / 2
            mid_idx = np.argmin(np.abs(vel_km2 - mid_km))
            vel_slice = vel_arr2[-n_2cycles:, mid_idx]
            displacement_m = np.cumsum(vel_slice * dt_sec)
            displacement_km = displacement_m / 1000.0

            ax_r3[3].bar(time_hours, vel_slice, width=time_hours[1] * 0.8 if len(time_hours) > 1 else 0.1,
                         color="steelblue", label="Velocity")
            ax_r3_twin = ax_r3[3].twinx()
            ax_r3_twin.plot(time_hours, displacement_km, "r--s", ms=3, label="Displacement")
            ax_r3_twin.set_ylabel("Displacement (km)")
            ax_r3[3].set(title=f"Transport at km {mid_km:.0f}",
                         xlabel="Time (h)", ylabel="Velocity (m/s)")
            lines1, labels1 = ax_r3[3].get_legend_handles_labels()
            lines2, labels2 = ax_r3_twin.get_legend_handles_labels()
            ax_r3[3].legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc="upper right")

            # R3C5: Eulerian vs Lagrangian salinity
            sal_at_fixed = sal_slice[:, mid_idx] if sal_arr.shape[1] > mid_idx else np.zeros(n_2cycles)
            particle_pos = np.clip(mid_km + displacement_km, sal_km[0], sal_km[-1])
            sal_lagrangian = [np.interp(p, sal_km, sal_slice[i])
                              for i, p in enumerate(particle_pos)]
            ax_r3[4].plot(time_hours, sal_at_fixed, "k-o", ms=2, label="Eulerian (fixed)")
            ax_r3[4].plot(time_hours, sal_lagrangian, "g->", ms=2, label="Lagrangian (parcel)")
            ax_r3[4].set(title=f"Salinity at km {mid_km:.0f}",
                         xlabel="Time (h)", ylabel="Salinity (PSU)")
            ax_r3[4].legend(fontsize=8)
            ax_r3[4].grid(True, alpha=0.3)
        else:
            _text_or_empty(ax_r3[3], "No Velocity Data")
            _text_or_empty(ax_r3[4], "No Velocity Data")
    else:
        for ax in ax_r3:
            _text_or_empty(ax, "No Salinity Data")

    # --- Save figure ---
    plt.subplots_adjust(left=0.04, right=0.98, top=0.95, bottom=0.05)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fig_path = fig_dir / "fig_hydrodynamics.png"
    fig.savefig(fig_path, dpi=250, bbox_inches="tight")
    fig.savefig(fig_dir / f"fig_hydrodynamics_{ts}.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # --- Metrics ---
    metrics_rows = []
    if obs_df is not None and depth_df is not None and "Tidal_range_m" in obs_df.columns:
        overall_mean, _, _ = _tidal_range_per_cycle(depth_df.values, None, tidal_period_steps)
        obs_vals = obs_df["Tidal_range_m"].dropna().values
        obs_km_v = obs_df.loc[obs_df["Tidal_range_m"].notna(), "Distance_km"].values
        sim_at_obs = np.interp(obs_km_v, depth_df.columns.astype(float).values, overall_mean)
        metrics_rows.append({
            "Variable": "Tidal_range", "Season": "all", "N": len(obs_vals),
            "RMSE": _rmse(obs_vals, sim_at_obs),
            "Pearson_R": _pearson_r(obs_vals, sim_at_obs),
            "Skill": _willmott_skill(obs_vals, sim_at_obs),
        })

    if obs_df is not None and sal_df is not None and "Salinity" in obs_df.columns:
        for sn in ("dry", "wet", "all"):
            s_mask = _get_season_masks(sal_df.index)
            if sn == "dry":
                mask = s_mask[0]
            elif sn == "wet":
                mask = s_mask[1]
            else:
                mask = np.ones(len(sal_df), dtype=bool)
            if not mask.any():
                continue
            sim_p = sal_df.loc[mask].mean(axis=0)
            obs_s = obs_df if sn == "all" else obs_df[obs_df.get("Season", pd.Series()) == sn]
            obs_s = obs_s.dropna(subset=["Salinity"])
            if len(obs_s) < 3:
                continue
            sim_at_obs = np.interp(obs_s["Distance_km"].values,
                                   sim_p.index.astype(float).values, sim_p.values)
            metrics_rows.append({
                "Variable": "Salinity", "Season": sn, "N": len(obs_s),
                "RMSE": _rmse(obs_s["Salinity"].values, sim_at_obs),
                "Pearson_R": _pearson_r(obs_s["Salinity"].values, sim_at_obs),
                "Skill": _willmott_skill(obs_s["Salinity"].values, sim_at_obs),
            })

    metrics_df = pd.DataFrame(metrics_rows)
    if not metrics_df.empty:
        metrics_df.to_csv(tbl_dir / "validation_hydrodynamics.csv",
                          index=False, float_format="%.4f")

    # Console
    print(f"\n{'='*65}")
    print("  HYDRODYNAMICS & TRANSPORT DIAGNOSTICS")
    print(f"{'='*65}")
    if metrics_df.empty:
        print("  Model-only (no observation data)")
    else:
        for _, row in metrics_df.iterrows():
            r_str = f"R={row['Pearson_R']:.3f}" if not np.isnan(row['Pearson_R']) else "R=---"
            s_str = f"Skill={row['Skill']:.3f}" if 'Skill' in row and not np.isnan(row.get('Skill', np.nan)) else ""
            print(f"  {row['Variable']:15s} ({row['Season']:4s})  "
                  f"N={int(row['N']):4d}  RMSE={row['RMSE']:7.3f}  {r_str}  {s_str}")
    print(f"  Figure → {fig_path}")
    print(f"{'='*65}")


if __name__ == "__main__":
    DEFAULT_OBS = Path("INPUT/Validation/obs_tidal_range.csv")
    parser = argparse.ArgumentParser(description="Hydrodynamics & transport diagnostics")
    parser.add_argument("--output-dir", default="OUT", help="Model output directory")
    parser.add_argument("--with-obs", action="store_true",
                        help=f"Enable obs overlay using default file ({DEFAULT_OBS})")
    parser.add_argument("--obs-file", default=None,
                        help="Observation CSV (Distance_km, Tidal_range_m[, Salinity, Season])")
    parser.add_argument("--save-dir", default=None, help="Output directory for figures/tables")
    args = parser.parse_args()
    obs = args.obs_file or (str(DEFAULT_OBS) if args.with_obs else None)
    plot(Path(args.output_dir), obs,
         args.save_dir or args.output_dir)
