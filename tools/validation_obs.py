#!/usr/bin/env python3
"""Shared observation helpers for validation plotting tools.

Supports both:
1. Compact model-unit observation tables used by the streamlined tools
2. Legacy `WaterQuality.csv` format used by the phase3/phase4 validation scripts
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_DRY_MONTHS = (11, 12, 1, 2, 3, 4)

VAR_SPECS = {
    "Sal": {"legacy_cols": ["Salinity"], "factor": 1.0},
    "SPM": {"legacy_cols": ["TSS (mg/L)"], "factor": 1.0},
    "O2": {"legacy_cols": ["DO (mg/L)"], "factor": 1000.0 / 32.0},
    "TOC": {"legacy_cols": ["TOC (mgC/L)"], "factor": 1000.0 / 12.01},
    "NH4": {"legacy_cols": ["NH4 (mgN/L)"], "factor": 1000.0 / 14.01},
    "NO3": {"legacy_cols": ["NO3 (mgN/L)"], "factor": 1000.0 / 14.01},
    "PO4": {"legacy_cols": ["PO4 (mgP/L)"], "factor": 1000.0 / 30.97},
    "Phy1": {"legacy_cols": ["Chl-a (ug/L)"], "factor": 15.0 / 12.01},
    "Si": {"legacy_cols": ["DSi (mgSi/L)"], "factor": 1000.0 / 28.09},
    "DIC": {"legacy_cols": ["DIC (mgC/L)"], "factor": 1000.0 / 12.01},
    "AT": {"legacy_cols": ["Alkalinity (mg/L CaCO3)"], "factor": 1000.0 / 50.04345},
    "PH": {"legacy_cols": ["pH"], "factor": 1.0},
    "pCO2": {"legacy_cols": ["pCO2 (ppmv)", "pCO2 (uatm)"], "factor": 1.0},
}

SIM_SCALES = {
    "SPM": 1000.0,
}

UNIT_LABEL_OVERRIDES = {
    "SPM": "g m⁻³",
}


def classify_validation_season(dates) -> np.ndarray:
    """Return dry/wet seasons using configurable tropical dry months.

    Override default months with environment variable `CGEM_DRY_MONTHS`,
    e.g. `CGEM_DRY_MONTHS=12,1,2,3,4,5`.
    """
    months = pd.DatetimeIndex(dates).month
    return np.where(np.isin(months, get_dry_months()), "dry", "wet")


def get_dry_months() -> tuple[int, ...]:
    """Return configured dry months for generic tropical-estuary season splitting."""
    raw = os.getenv("CGEM_DRY_MONTHS", "").strip()
    if not raw:
        return DEFAULT_DRY_MONTHS
    try:
        months = tuple(int(part.strip()) for part in raw.split(",") if part.strip())
    except ValueError:
        return DEFAULT_DRY_MONTHS
    valid = tuple(month for month in months if 1 <= month <= 12)
    return valid if valid else DEFAULT_DRY_MONTHS


def normalize_season_labels(series: pd.Series | None) -> pd.Series:
    """Normalize mixed season labels to 'dry'/'wet'."""
    if series is None:
        return pd.Series(dtype=object)
    raw = series.astype(str).str.strip().str.lower()
    return pd.Series(
        np.where(raw.str.contains("dry"), "dry",
                 np.where(raw.str.contains("wet|rain"), "wet", "unknown")),
        index=series.index,
    )


def default_obs_path(kind: str) -> Path:
    """Return the preferred default observation file for a tool kind."""
    candidates = [
        Path("INPUT/Validation/obs_longitudinal_profiles.csv") if kind == "profiles"
        else Path("INPUT/Validation/obs_timeseries.csv"),
        Path("INPUT/Validation/WaterQuality.csv"),
        Path("OUT/PlotINPUT/WaterQuality.csv"),
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[-1]


def get_sim_scale(var_name: str) -> float:
    return float(SIM_SCALES.get(var_name, 1.0))


def get_unit_override(var_name: str) -> str | None:
    return UNIT_LABEL_OVERRIDES.get(var_name)


def _ensure_base_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rename_map = {}
    if "Location" in out.columns and "Distance_km" not in out.columns:
        rename_map["Location"] = "Distance_km"
    if "Site" in out.columns and "Station" not in out.columns:
        rename_map["Site"] = "Station"
    if rename_map:
        out = out.rename(columns=rename_map)

    if "Date" in out.columns:
        out["Date"] = pd.to_datetime(out["Date"], errors="coerce")
    if "Distance_km" in out.columns:
        out["Distance_km"] = pd.to_numeric(out["Distance_km"], errors="coerce")
    return out


def _is_compact_format(df: pd.DataFrame) -> bool:
    return "Distance_km" in df.columns and any(var in df.columns for var in VAR_SPECS)


def _standardize_compact(df: pd.DataFrame) -> pd.DataFrame:
    out = _ensure_base_columns(df)
    if "Season" in out.columns:
        out["Season"] = normalize_season_labels(out["Season"])
    elif "Date" in out.columns:
        out["Season"] = normalize_season_labels(
            pd.Series(classify_validation_season(out["Date"].dropna()), index=out["Date"].dropna().index)
            .reindex(out.index)
        )
    return out


def _standardize_legacy(df: pd.DataFrame) -> pd.DataFrame:
    out = _ensure_base_columns(df)
    keep_cols = [col for col in ("Date", "Station", "Distance_km", "Source") if col in out.columns]
    standardized = out[keep_cols].copy() if keep_cols else pd.DataFrame(index=out.index)

    if "Season" in out.columns:
        standardized["Season"] = normalize_season_labels(out["Season"])
    elif "Date" in out.columns:
        season_labels = pd.Series(classify_validation_season(out["Date"].dropna()), index=out["Date"].dropna().index)
        standardized["Season"] = normalize_season_labels(season_labels.reindex(out.index))

    for var_name, spec in VAR_SPECS.items():
        for legacy_col in spec["legacy_cols"]:
            if legacy_col in out.columns:
                standardized[var_name] = pd.to_numeric(out[legacy_col], errors="coerce") * float(spec["factor"])
                break

    return standardized


def standardize_observations(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize observation tables to a common tool-friendly schema."""
    base = _ensure_base_columns(df)
    if _is_compact_format(base):
        return _standardize_compact(base)
    return _standardize_legacy(base)


def load_standardized_observations(obs_file: str | None = None, kind: str = "profiles") -> tuple[pd.DataFrame | None, Path]:
    """Load observations and normalize them to the shared schema."""
    target = Path(obs_file) if obs_file else default_obs_path(kind)
    if not target.exists():
        return None, target
    df = pd.read_csv(target)
    return standardize_observations(df), target