#!/usr/bin/env python3
"""
validate_saigon.py — Comprehensive C-GEM validation for the Saigon River.

Generates:
  1. Longitudinal profile validation (model seasonal-mean vs obs)
  2. Salinity & tidal range longitudinal validation
  3. Seasonal time-series validation at monitoring stations
  4. Carbonate system longitudinal validation (pCO2, DIC, pH, AT)
  5. Validation metrics table (RMSE, Bias, MAE, R², n)

All figures → OUT/figures/
All tables  → OUT/tables/

Usage:
    python validate_saigon.py [--out-dir OUT]
"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "OUT"
VAL_DIR = ROOT / "INPUT" / "Validation"
OBS_PROFILES = VAL_DIR / "obs_longitudinal_profiles.csv"
OBS_TIMESERIES = VAL_DIR / "obs_timeseries.csv"
OBS_TIDAL = VAL_DIR / "obs_tidal_range.csv"
GEOM_FILE = OUT_DIR / "geometry.csv"
FIG_DIR = OUT_DIR / "figures"
TBL_DIR = OUT_DIR / "tables"

# ── Constants ──────────────────────────────────────────────────────────
MW_N  = 14.007    # g/mol
MW_P  = 30.974
MW_C  = 12.011
MW_O2 = 32.0
MW_Si = 28.085

# C:Chla mass ratio (default 50)
CCHLA_RATIO = 50.0
CHLA_CONV = MW_C / CCHLA_RATIO   # mmol C/m³ → µg Chla/L

DRY_MONTHS  = {11, 12, 1, 2, 3, 4}
WET_MONTHS  = {5, 6, 7, 8, 9, 10}
SEASON_COLORS = {"dry": "#D84315", "wet": "#1565C0"}

# Station mapping from WaterQuality.csv
STATIONS_WQ = {
    "Nga Bay": 0, "Vam Sat": 24, "Nha Be": 44,
    "Bach Dang": 70, "Binh Loi": 88, "Phu Cuong": 114, "Ben Suc": 156,
}



# ── Variable definitions ──────────────────────────────────────────────
# Key: model var name
# Value: (obs_column_in_WQ, display_label, unit_label, conv_factor)
WQ_VARS = {
    "Sal":  ("Salinity_PSU",   "Salinity",   "PSU",      1.0),
    "Chla": ("Chla_ugL",      "Chl-a",      "µg/L",     None),  # special
    "O2":   ("DO_mgL",        "DO",         "mg O₂/L",  MW_O2 / 1000),
    "NH4":  ("NH4_mgNL",      "NH₄",        "mg N/L",   MW_N  / 1000),
    "NO3":  ("NO3_mgNL",      "NO₃",        "mg N/L",   MW_N  / 1000),
    "PO4":  ("PO4_mgPL",      "PO₄",        "mg P/L",   MW_P  / 1000),
    "TOC":  ("TOC_mgCL",      "TOC",        "mg C/L",   MW_C  / 1000),
    "SPM":  ("TSS_mgL",       "TSS",        "mg/L",     1000.0),  # model: kg/m³ → mg/L
    "Si":   ("DSi_mgSiL",     "DSi",        "mg Si/L",  MW_Si / 1000),
}

CARB_VARS = {
    "pCO2": ("pCO2_uatm",            "pCO₂",  "µatm",        1.0),
    "DIC":  ("DIC_mgCL",             "DIC",   "mg C/L",      MW_C / 1000),
    "PH":   ("pH",                   "pH",    "—",           1.0),
    "AT":   ("Alkalinity_mgCaCO3L",  "AT",    "mg CaCO₃/L",  0.050045),
    # AT model: µeq/L → mg CaCO₃/L: CaCO₃ divalent → 1 meq = 50.045 mg CaCO₃
}

# ── Helper functions ──────────────────────────────────────────────────
def classify_season(dt):
    return "dry" if dt.month in DRY_MONTHS else "wet"


def load_geometry():
    df = pd.read_csv(GEOM_FILE)
    df["km"] = df["Distance"] / 1000.0
    return df


def load_model(var_name):
    """Load model output CSV. Returns DataFrame with datetime index."""
    path = OUT_DIR / f"{var_name}.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path, parse_dates=["Date"], index_col="Date")
    return df


def load_model_chla():
    """Combine Phy1+Phy2 and convert to Chla (µg/L)."""
    phy1 = load_model("Phy1")
    phy2 = load_model("Phy2")
    if phy1 is None or phy2 is None:
        return None
    return phy1.add(phy2, fill_value=0) * CHLA_CONV


def load_obs_profiles():
    """Load longitudinal profile observations (WQ + carbonate)."""
    df = pd.read_csv(OBS_PROFILES)
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    df["Distance_km"] = pd.to_numeric(df["Distance_km"], errors="coerce")
    df = df.dropna(subset=["Date", "Distance_km"])
    df["season"] = df["Season"].str.strip().str.lower().apply(
        lambda x: "dry" if "dry" in str(x) else "wet"
    )
    return df


def load_obs_timeseries():
    """Load time-series observations at fixed stations."""
    if not OBS_TIMESERIES.exists():
        return None
    df = pd.read_csv(OBS_TIMESERIES)
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    df["Distance_km"] = pd.to_numeric(df["Distance_km"], errors="coerce")
    df = df.dropna(subset=["Date", "Distance_km"])
    df["season"] = df["Season"].str.strip().str.lower().apply(
        lambda x: "dry" if "dry" in str(x) else "wet"
    )
    return df


def seasonal_profile(model_df, geom, season, conv=1.0):
    """Compute seasonal-mean, p5, p95 longitudinal profile."""
    seasons = model_df.index.map(classify_season)
    s_df = model_df.loc[seasons == season]
    if s_df.empty:
        n = len(geom)
        nan = np.full(n, np.nan)
        return geom["km"].values, nan, nan, nan
    km = geom["km"].values
    vals = s_df.values * conv
    n = min(len(km), vals.shape[1])
    mean_v = np.nanmean(vals[:, :n], axis=0)
    p5  = np.nanpercentile(vals[:, :n], 5, axis=0)
    p95 = np.nanpercentile(vals[:, :n], 95, axis=0)
    return km[:n], mean_v, p5, p95


def find_cell(geom, km_target):
    """Return (cell_index_0based, cell_col_name) for nearest cell."""
    idx = np.argmin(np.abs(geom["km"].values - km_target))
    return idx, f"Cell_{idx + 1}"


def _shade_urban(ax, km_range=(40, 110)):
    """Shade urban corridor on longitudinal plots."""
    ax.axvspan(*km_range, alpha=0.06, color="#9C27B0", zorder=0)


def _common_legend(fig, n_cols=4, loc="lower center", bbox=(0.5, 0.002)):
    handles = [
        Line2D([0], [0], color=SEASON_COLORS["dry"], lw=2, label="Model (dry)"),
        Line2D([0], [0], color=SEASON_COLORS["wet"], lw=2, label="Model (wet)"),
        Line2D([0], [0], marker="o", ls="", color=SEASON_COLORS["dry"],
               ms=5, label="Obs (dry)"),
        Line2D([0], [0], marker="o", ls="", color=SEASON_COLORS["wet"],
               ms=5, label="Obs (wet)"),
    ]
    fig.legend(handles=handles, loc=loc, ncol=n_cols, fontsize=9,
               frameon=True, bbox_to_anchor=bbox)


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 1: Longitudinal WQ Profiles
# ═══════════════════════════════════════════════════════════════════════
def fig_longitudinal_profiles():
    print("[Fig 1] Longitudinal WQ profiles...")
    geom = load_geometry()
    obs = load_obs_profiles()

    vars_to_plot = ["Sal", "Chla", "O2", "NH4", "NO3", "PO4", "TOC", "SPM", "Si"]
    n = len(vars_to_plot)
    fig, axes = plt.subplots(3, 3, figsize=(16, 13))
    axes = axes.flatten()

    for i, var in enumerate(vars_to_plot):
        ax = axes[i]
        _shade_urban(ax)

        obs_col, label, unit, conv = WQ_VARS[var]

        # Model profiles
        for season in ("dry", "wet"):
            if var == "Chla":
                chla_df = load_model_chla()
                if chla_df is not None:
                    km, mean_v, p5, p95 = seasonal_profile(chla_df, geom, season, conv=1.0)
                else:
                    continue
            else:
                m_df = load_model(var)
                if m_df is None:
                    continue
                km, mean_v, p5, p95 = seasonal_profile(m_df, geom, season, conv=conv)

            ax.plot(km, mean_v, color=SEASON_COLORS[season], lw=1.5,
                    label=f"Model ({season})")
            ax.fill_between(km, p5, p95, color=SEASON_COLORS[season], alpha=0.10)

        # Observations
        if obs_col in obs.columns:
            for season in ("dry", "wet"):
                s = obs[obs["season"] == season].copy()
                s[obs_col] = pd.to_numeric(s[obs_col], errors="coerce")
                s = s.dropna(subset=[obs_col])
                ax.scatter(s["Distance_km"], s[obs_col], s=15,
                           color=SEASON_COLORS[season], edgecolor="none",
                           alpha=0.55, zorder=5)

        ax.set_ylabel(f"{label} ({unit})", fontsize=10)
        ax.set_xlabel("Distance from mouth (km)", fontsize=9)
        ax.set_xlim(0, 210)
        ax.grid(True, alpha=0.15)
        ax.set_title(label, fontweight="bold", fontsize=11)

    _common_legend(fig)
    fig.suptitle("Longitudinal Profile Validation — Water Quality",
                 fontweight="bold", fontsize=14, y=0.99)
    plt.tight_layout(rect=(0, 0.035, 1, 0.97))
    path = FIG_DIR / "fig_validation_longitudinal_WQ.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {path}")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 2: Salinity & Tidal Range
# ═══════════════════════════════════════════════════════════════════════
def fig_salinity_tidal():
    print("[Fig 2] Salinity & tidal range...")
    geom = load_geometry()
    obs_profiles = load_obs_profiles()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # --- Panel A: Salinity longitudinal ---
    _shade_urban(ax1)
    sal_df = load_model("Sal")
    if sal_df is not None:
        for season in ("dry", "wet"):
            km, mean_v, p5, p95 = seasonal_profile(sal_df, geom, season)
            ax1.plot(km, mean_v, color=SEASON_COLORS[season], lw=1.5,
                     label=f"Model ({season})")
            ax1.fill_between(km, p5, p95, color=SEASON_COLORS[season], alpha=0.10)

    if "Salinity_PSU" in obs_profiles.columns:
        for season in ("dry", "wet"):
            s = obs_profiles[obs_profiles["season"] == season].copy()
            s["Salinity_PSU"] = pd.to_numeric(s["Salinity_PSU"], errors="coerce")
            s = s.dropna(subset=["Salinity_PSU"])
            ax1.scatter(s["Distance_km"], s["Salinity_PSU"], s=20,
                        color=SEASON_COLORS[season], edgecolor="none",
                        alpha=0.6, zorder=5, label=f"Obs ({season})")

    ax1.set_xlabel("Distance from mouth (km)")
    ax1.set_ylabel("Salinity (PSU)")
    ax1.set_xlim(0, 210)
    ax1.set_title("(a) Salinity", fontweight="bold")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.15)

    # --- Panel B: Tidal range ---
    _shade_urban(ax2)
    wd_df = load_model("waterDepth")
    if wd_df is not None:
        # Model tidal range: daily max - min water depth per cell
        daily_max = wd_df.resample("D").max()
        daily_min = wd_df.resample("D").min()
        tidal_range = daily_max - daily_min
        km = geom["km"].values
        mean_tr = tidal_range.mean(axis=0).values
        n = min(len(km), len(mean_tr))
        ax2.plot(km[:n], mean_tr[:n], "k-", lw=1.5, label="Model (mean)")

    if OBS_TIDAL.exists():
        tidal_obs = pd.read_csv(OBS_TIDAL)
        tidal_obs["Distance_km"] = pd.to_numeric(tidal_obs["Distance_km"], errors="coerce")
        tidal_obs["Tidal_range_m"] = pd.to_numeric(tidal_obs["Tidal_range_m"],
                                                      errors="coerce")
        tidal_obs = tidal_obs.dropna()
        # Box plot at each station
        locs = sorted(tidal_obs["Distance_km"].unique())
        for loc in locs:
            subset = tidal_obs[tidal_obs["Distance_km"] == loc]["Tidal_range_m"]
            bp = ax2.boxplot(subset, positions=[loc], widths=6,
                             patch_artist=True, zorder=5,
                             boxprops=dict(facecolor="#90CAF9", alpha=0.7),
                             medianprops=dict(color="#D84315", lw=1.5),
                             flierprops=dict(marker=".", ms=3))

    ax2.set_xlabel("Distance from mouth (km)")
    ax2.set_ylabel("Tidal range (m)")
    ax2.set_xlim(0, 210)
    ax2.set_title("(b) Tidal Range", fontweight="bold")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.15)

    fig.suptitle("Salinity & Tidal Range Validation", fontweight="bold", fontsize=13)
    plt.tight_layout()
    path = FIG_DIR / "fig_validation_salinity_tidal.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {path}")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 3: Seasonal Time Series at Monitoring Stations
# ═══════════════════════════════════════════════════════════════════════
def fig_seasonal_timeseries():
    print("[Fig 3] Seasonal time series at monitoring stations...")
    geom = load_geometry()
    obs_ts = load_obs_timeseries()

    # Focus stations (4 key locations)
    stations = [
        ("Nha Be (44 km)", 44),
        ("Bach Dang (70 km)", 70),
        ("Phu Cuong (114 km)", 114),
        ("Ben Suc (156 km)", 156),
    ]
    vars_ts = ["Sal", "Chla", "O2", "NH4", "NO3", "PO4"]

    fig, axes = plt.subplots(len(vars_ts), len(stations), figsize=(18, 20),
                              squeeze=False)

    for col, (stn_name, stn_km) in enumerate(stations):
        cell_idx, cell_col = find_cell(geom, stn_km)

        for row, var in enumerate(vars_ts):
            ax = axes[row][col]

            obs_col, label, unit, conv = WQ_VARS[var]

            # Model time series (monthly mean)
            if var == "Chla":
                m_df = load_model_chla()
                conv_eff = 1.0
            else:
                m_df = load_model(var)
                conv_eff = conv

            if m_df is not None and cell_col in m_df.columns:
                ts = m_df[cell_col] * conv_eff
                monthly = ts.resample("MS").mean()
                ax.plot(monthly.index, monthly.values, "k-", lw=1.0,
                        label="Model", zorder=3)

                # Shade wet season bands
                for yr in sorted(monthly.index.year.unique()):
                    ws = pd.Timestamp(f"{yr}-05-01")
                    we = pd.Timestamp(f"{yr}-10-31")
                    ax.axvspan(ws, we, alpha=0.05, color="#1565C0", zorder=0)

            # Observations at this station (±8 km)
            if obs_ts is not None and obs_col in obs_ts.columns:
                stn_obs = obs_ts[
                    (obs_ts["Distance_km"] >= stn_km - 8) &
                    (obs_ts["Distance_km"] <= stn_km + 8)
                ].copy()
                stn_obs[obs_col] = pd.to_numeric(stn_obs[obs_col], errors="coerce")
                stn_obs = stn_obs.dropna(subset=[obs_col])

                # Split into model-period and extended
                model_start = m_df.index.min() if m_df is not None else pd.Timestamp("2017-01-01")
                model_end = m_df.index.max() if m_df is not None else pd.Timestamp("2018-12-31")

                in_period = stn_obs[
                    (stn_obs["Date"] >= model_start) &
                    (stn_obs["Date"] <= model_end)
                ]
                ext_period = stn_obs[
                    (stn_obs["Date"] < model_start) |
                    (stn_obs["Date"] > model_end)
                ]

                for season in ("dry", "wet"):
                    sd = in_period[in_period["season"] == season]
                    if not sd.empty:
                        ax.scatter(sd["Date"], sd[obs_col], s=15,
                                   color=SEASON_COLORS[season], edgecolor="none",
                                   alpha=0.7, zorder=5)

                if not ext_period.empty:
                    ax.scatter(ext_period["Date"], ext_period[obs_col],
                               s=12, marker="^", color="#757575",
                               edgecolor="none", alpha=0.4, zorder=4)

            ax.set_ylabel(f"{label} ({unit})", fontsize=8)
            ax.grid(True, alpha=0.15)
            ax.tick_params(labelsize=7)

            if row == 0:
                ax.set_title(stn_name, fontweight="bold", fontsize=10)
            if row < len(vars_ts) - 1:
                ax.set_xticklabels([])
            else:
                ax.xaxis.set_major_locator(mdates.YearLocator())
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
                ax.tick_params(axis="x", rotation=30, labelsize=7)

    # Legend
    handles = [
        Line2D([0], [0], color="k", lw=1.2, label="Model (monthly)"),
        Line2D([0], [0], marker="o", ls="", color=SEASON_COLORS["dry"],
               ms=5, label="Obs-dry"),
        Line2D([0], [0], marker="o", ls="", color=SEASON_COLORS["wet"],
               ms=5, label="Obs-wet"),
        Line2D([0], [0], marker="^", ls="", color="#757575",
               ms=5, label="SDNR extended"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9,
               frameon=True, bbox_to_anchor=(0.5, 0.002))
    fig.suptitle("Seasonal Time-Series Validation", fontweight="bold",
                 fontsize=14, y=0.995)
    plt.tight_layout(rect=(0, 0.03, 1, 0.98))
    path = FIG_DIR / "fig_validation_seasonal_timeseries.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {path}")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 4: Carbonate Longitudinal Profiles
# ═══════════════════════════════════════════════════════════════════════
def fig_carbonate_profiles():
    print("[Fig 4] Carbonate longitudinal profiles...")
    geom = load_geometry()
    obs_profiles = load_obs_profiles()
    carb_has_data = obs_profiles is not None and any(
        c in obs_profiles.columns for c in ["pCO2_uatm", "DIC_mgCL", "pH", "Alkalinity_mgCaCO3L"]
    )
    if not carb_has_data:
        print("  SKIP: no carbonate data in profile observations")
        return

    carb_vars = ["pCO2", "DIC", "PH", "AT"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for i, var in enumerate(carb_vars):
        ax = axes[i]
        _shade_urban(ax)
        obs_col, label, unit, conv = CARB_VARS[var]

        # Fix AT conversion: model AT is in mmol/m³ (= µeq/L)
        # 1 mmol CaCO₃ = 100.09 mg → conv = 100.09/1000 for mg CaCO₃/L
        conv_eff = conv

        m_df = load_model(var)
        if m_df is not None:
            for season in ("dry", "wet"):
                km, mean_v, p5, p95 = seasonal_profile(m_df, geom, season, conv=conv_eff)
                ax.plot(km, mean_v, color=SEASON_COLORS[season], lw=1.5,
                        label=f"Model ({season})")
                ax.fill_between(km, p5, p95, color=SEASON_COLORS[season], alpha=0.10)

        if obs_col in obs_profiles.columns:
            for season in ("dry", "wet"):
                s = obs_profiles[obs_profiles["season"] == season].copy()
                s[obs_col] = pd.to_numeric(s[obs_col], errors="coerce")
                s = s.dropna(subset=[obs_col])
                ax.scatter(s["Distance_km"], s[obs_col], s=25,
                           color=SEASON_COLORS[season], edgecolor="none",
                           alpha=0.6, zorder=5, marker="s")

        ax.set_xlabel("Distance from mouth (km)", fontsize=9)
        ax.set_ylabel(f"{label} ({unit})", fontsize=10)
        ax.set_xlim(0, 210)
        ax.set_title(label, fontweight="bold", fontsize=11)
        ax.grid(True, alpha=0.15)
        if i == 0:
            ax.legend(fontsize=8)

    fig.suptitle("Carbonate System Longitudinal Validation",
                 fontweight="bold", fontsize=13)
    _common_legend(fig)
    plt.tight_layout(rect=(0, 0.04, 1, 0.96))
    path = FIG_DIR / "fig_validation_carbonate_longitudinal.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {path}")


# ═══════════════════════════════════════════════════════════════════════
# TABLE: Validation Metrics
# ═══════════════════════════════════════════════════════════════════════
def compute_metrics(obs_vals, mod_vals):
    """Compute RMSE, Bias, MAE, R² between matched arrays."""
    mask = np.isfinite(obs_vals) & np.isfinite(mod_vals)
    obs_v = obs_vals[mask]
    mod_v = mod_vals[mask]
    n = len(obs_v)
    if n < 2:
        return {"n": n, "RMSE": np.nan, "Bias": np.nan, "MAE": np.nan, "R2": np.nan}
    diff = mod_v - obs_v
    rmse = np.sqrt(np.mean(diff ** 2))
    bias = np.mean(diff)
    mae  = np.mean(np.abs(diff))
    ss_res = np.sum(diff ** 2)
    ss_tot = np.sum((obs_v - np.mean(obs_v)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return {"n": n, "RMSE": rmse, "Bias": bias, "MAE": mae, "R2": r2}


def table_validation_metrics():
    print("[Table] Validation metrics...")
    geom = load_geometry()
    obs_profiles = load_obs_profiles()

    all_stations = list(STATIONS_WQ.items())
    rows = []

    # WQ variables
    for var, (obs_col, label, unit, conv) in WQ_VARS.items():
        if obs_col not in obs_profiles.columns:
            continue

        if var == "Chla":
            m_df = load_model_chla()
            conv_eff = 1.0
        else:
            m_df = load_model(var)
            conv_eff = conv

        if m_df is None:
            continue

        for stn_name, stn_km in all_stations:
            _, cell_col = find_cell(geom, stn_km)
            stn_obs = obs_profiles[
                (obs_profiles["Distance_km"] >= stn_km - 5) &
                (obs_profiles["Distance_km"] <= stn_km + 5)
            ].copy()
            stn_obs[obs_col] = pd.to_numeric(stn_obs[obs_col], errors="coerce")
            stn_obs = stn_obs.dropna(subset=[obs_col])

            if len(stn_obs) < 3:
                continue

            obs_vals = []
            mod_vals = []
            for _, r in stn_obs.iterrows():
                cidx = m_df.index.get_indexer([r["Date"]], method="nearest")[0]
                if 0 <= cidx < len(m_df):
                    mod_vals.append(m_df.iloc[cidx][cell_col] * conv_eff)
                    obs_vals.append(r[obs_col])

            obs_arr = np.array(obs_vals, dtype=float)
            mod_arr = np.array(mod_vals, dtype=float)
            met = compute_metrics(obs_arr, mod_arr)

            for season in ("all", "dry", "wet"):
                if season == "all":
                    o, m = obs_arr, mod_arr
                else:
                    s_mask = np.array([classify_season(d) == season for d in stn_obs["Date"]])
                    if len(s_mask) != len(obs_arr):
                        s_mask = s_mask[:len(obs_arr)]
                    o = obs_arr[s_mask]
                    m = mod_arr[s_mask]
                met = compute_metrics(o, m)
                rows.append({
                    "Variable": label, "Unit": unit,
                    "Station": stn_name, "Distance_km": stn_km,
                    "Season": season, **met
                })

    # Carbonate variables (also in obs_longitudinal_profiles)
    for var, (obs_col, label, unit, conv) in CARB_VARS.items():
        if obs_col not in obs_profiles.columns:
            continue
        m_df = load_model(var)
        if m_df is None:
            continue
        conv_eff = conv

        for stn_name, stn_km in all_stations:
            _, cell_col = find_cell(geom, stn_km)
            stn_obs = obs_profiles[
                (obs_profiles["Distance_km"] >= stn_km - 8) &
                (obs_profiles["Distance_km"] <= stn_km + 8)
            ].copy()
            stn_obs[obs_col] = pd.to_numeric(stn_obs[obs_col], errors="coerce")
            stn_obs = stn_obs.dropna(subset=[obs_col])

            if len(stn_obs) < 2:
                continue

            # For carbonate: use model seasonal mean at nearest cell
            obs_vals = stn_obs[obs_col].values.astype(float)
            mod_vals = []
            for _, r in stn_obs.iterrows():
                cidx = m_df.index.get_indexer([r["Date"]], method="nearest")[0]
                if 0 <= cidx < len(m_df):
                    mod_vals.append(m_df.iloc[cidx][cell_col] * conv_eff)
                else:
                    mod_vals.append(np.nan)
            mod_arr = np.array(mod_vals, dtype=float)
            obs_arr = obs_vals[:len(mod_arr)]

            for season in ("all", "dry", "wet"):
                if season == "all":
                    o, m = obs_arr, mod_arr
                else:
                    s_mask = np.array([
                        classify_season(d) == season for d in stn_obs["Date"]
                    ])
                    s_mask = s_mask[:len(obs_arr)]
                    o = obs_arr[s_mask]
                    m = mod_arr[s_mask]
                met = compute_metrics(o, m)
                if met["n"] >= 2:
                    rows.append({
                        "Variable": label, "Unit": unit,
                        "Station": stn_name, "Distance_km": stn_km,
                        "Season": season, **met
                    })

    if not rows:
        print("  No metrics computed")
        return

    df = pd.DataFrame(rows)
    df = df.sort_values(["Variable", "Distance_km", "Season"])

    # Full table
    full_path = TBL_DIR / "validation_metrics_full.csv"
    df.to_csv(full_path, index=False, float_format="%.4f")
    print(f"  → {full_path}  ({len(df)} rows)")

    # Summary table (all seasons, aggregated)
    summary = df[df["Season"] == "all"].copy()
    summary_path = TBL_DIR / "validation_metrics_summary.csv"
    summary.to_csv(summary_path, index=False, float_format="%.4f")
    print(f"  → {summary_path}  ({len(summary)} rows)")

    # Print summary
    print("\n  ┌─────────────────────────────────────────────────────────────────┐")
    print("  │  VALIDATION SUMMARY (all seasons)                              │")
    print("  ├─────────┬─────────────┬──────┬────────┬─────────┬──────────────┤")
    print("  │ Var     │ Station     │  n   │  RMSE  │  Bias   │   R²         │")
    print("  ├─────────┼─────────────┼──────┼────────┼─────────┼──────────────┤")
    for _, r in summary.iterrows():
        print(f"  │ {r['Variable']:7s} │ {r['Station']:11s} │ {r['n']:4.0f} │"
              f" {r['RMSE']:6.3f} │ {r['Bias']:+7.3f} │ {r['R2']:+7.3f}      │")
    print("  └─────────┴─────────────┴──────┴────────┴─────────┴──────────────┘")

    return df


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════
def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TBL_DIR.mkdir(parents=True, exist_ok=True)

    if not GEOM_FILE.exists():
        print("ERROR: Run the model first — no geometry.csv in OUT/")
        return 1
    if not OBS_PROFILES.exists():
        print(f"ERROR: Observation file missing: {OBS_PROFILES}")
        return 1

    print(f"Root:     {ROOT}")
    print(f"Outputs:  {OUT_DIR}")
    print(f"Figures → {FIG_DIR}")
    print(f"Tables  → {TBL_DIR}")
    print()

    fig_longitudinal_profiles()
    fig_salinity_tidal()
    fig_seasonal_timeseries()
    fig_carbonate_profiles()
    table_validation_metrics()

    print("\n✓ All validation outputs generated.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
