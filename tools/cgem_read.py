"""
cgem_read.py — I/O and utility module for C-GEM model analysis.

Core capabilities:
  - Read model output (CSV / CGEMBIN binary) as DataFrames
  - Season classification (tropical dry/wet)
  - Variable metadata (names, units)
  - Geometry and parameter loading

Usage:
    from cgem_read import load_output, load_geometry, list_variables
    from cgem_read import classify_season, get_variable_info
"""

import struct
import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MAGIC = b"CGEMBIN"
_HEADER_V2 = struct.Struct("<8s 3I 2d")  # magic, version, cells, vals/rec, dt, dx

# Default output directory (relative to project root)
DEFAULT_OUT_DIR = Path("OUT")

# Variable metadata: name -> (long_name, unit)
VARIABLE_INFO = {
    "Phy1":     ("Diatoms",                        "mmol C m⁻³"),
    "Phy2":     ("Non-siliceous phytoplankton",     "mmol C m⁻³"),
    "Si":       ("Dissolved silica",                "mmol Si m⁻³"),
    "NO3":      ("Nitrate",                         "mmol N m⁻³"),
    "NH4":      ("Ammonium",                        "mmol N m⁻³"),
    "PO4":      ("Phosphate",                       "mmol P m⁻³"),
    "PIP":      ("Particulate inorganic phosphorus","mmol P m⁻³"),
    "O2":       ("Dissolved oxygen",                "mmol O₂ m⁻³"),
    "TOC":      ("Total organic carbon",            "mmol C m⁻³"),
    "Sal":      ("Salinity",                        "PSU"),
    "SPM":      ("Suspended particulate matter",    "kg m⁻³"),
    "DIC":      ("Dissolved inorganic carbon",      "mmol C m⁻³"),
    "AT":       ("Total alkalinity",                "µeq L⁻¹"),
    "pCO2":     ("Partial pressure of CO₂",         "µatm"),
    "PH":       ("pH",                              "—"),
    "CO2":      ("Dissolved CO₂",                   "mmol C m⁻³"),
    # Hydrodynamics
    "velocity":   ("Flow velocity",                 "m s⁻¹"),
    "waterDepth": ("Water depth",                   "m"),
    "Discharge":  ("Discharge",                     "m³ s⁻¹"),
    "freeArea":   ("Free cross-section area",       "m²"),
    "totalArea":  ("Total cross-section area",      "m²"),
    "disp":       ("Dispersion coefficient",        "m² s⁻¹"),
    "tau_b":      ("Bed shear stress",              "Pa"),
    # Reaction rates
    "Reaction_NPP_NO3":          ("Net primary production (NO₃ uptake)",   "mmol m⁻³ s⁻¹"),
    "Reaction_NPP_NH4":          ("Net primary production (NH₄ uptake)",   "mmol m⁻³ s⁻¹"),
    "Reaction_phydeath":         ("Phytoplankton mortality",               "mmol m⁻³ s⁻¹"),
    "Reaction_adegradation":     ("Aerobic degradation",                   "mmol m⁻³ s⁻¹"),
    "Reaction_denitrification":  ("Denitrification",                       "mmol m⁻³ s⁻¹"),
    "Reaction_nitrification":    ("Nitrification",                         "mmol m⁻³ s⁻¹"),
    "Reaction_O2_exchange":      ("O₂ air-water exchange",                 "mmol m⁻³ s⁻¹"),
    "Reaction_CO2_exchange":     ("CO₂ air-water exchange",                "mmol m⁻³ s⁻¹"),
    "Reaction_DIC":              ("DIC source/sink",                       "mmol m⁻³ s⁻¹"),
    "Reaction_TA":               ("Alkalinity source/sink",                "mmol m⁻³ s⁻¹"),
    # Diagnostic limitation factors
    "Diag_fN":       ("Nitrogen limitation",     "—"),
    "Diag_fP":       ("Phosphorus limitation",   "—"),
    "Diag_fSi":      ("Silicon limitation",      "—"),
    "Diag_fI":       ("Light limitation",        "—"),
    "Diag_KD":       ("Light attenuation (Kd)",  "m⁻¹"),
}

# ---------------------------------------------------------------------------
# Season classification
# ---------------------------------------------------------------------------
# Season classification (tropical monsoon convention)
DRY_MONTHS = (11, 12, 1, 2, 3, 4)
WET_MONTHS = (5, 6, 7, 8, 9, 10)


def classify_season(dates) -> np.ndarray:
    """Classify dates as 'dry' or 'wet' (tropical monsoon convention).

    Parameters
    ----------
    dates : datetime-like, Series, or DatetimeIndex

    Returns
    -------
    numpy array of strings ('dry' or 'wet')
    """
    months = pd.DatetimeIndex(dates).month
    return np.where(np.isin(months, DRY_MONTHS), "dry", "wet")


# ---------------------------------------------------------------------------
# Binary reader
# ---------------------------------------------------------------------------

def _read_binary(path: Path) -> tuple:
    """Read a CGEMBIN file. Returns (cell_count, dt, dx, time_arr, data_2d)."""
    with open(path, "rb") as fh:
        raw = fh.read(_HEADER_V2.size)
        if len(raw) < _HEADER_V2.size:
            raise ValueError(f"{path}: too short for CGEMBIN header")

        magic, version, cell_count, vals_per_rec, dt, dx = _HEADER_V2.unpack(raw)
        if not magic.rstrip(b"\x00").startswith(MAGIC):
            raise ValueError(f"{path}: invalid CGEMBIN magic")

        # Handle v1 padding (4 extra bytes)
        header_size = _HEADER_V2.size
        if version <= 1:
            fh.read(4)  # skip padding
            header_size += 4

        record_size = vals_per_rec * 8
        file_size = os.fstat(fh.fileno()).st_size
        data_bytes = file_size - header_size
        n_records = data_bytes // record_size

        if n_records == 0:
            return cell_count, dt, dx, np.empty(0), np.empty((0, cell_count))

        buf = fh.read(record_size * n_records)
        data = np.frombuffer(buf, dtype="<f8").reshape(n_records, vals_per_rec)
        time_arr = data[:, 0].copy()
        values = data[:, 1:].copy()
        return cell_count, dt, dx, time_arr, values


# ---------------------------------------------------------------------------
# CSV reader
# ---------------------------------------------------------------------------

def _read_csv(path: Path) -> pd.DataFrame:
    """Read a C-GEM CSV output file. Returns DataFrame with Date index."""
    df = pd.read_csv(path, parse_dates=["Date"], index_col="Date")
    return df


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _detect_start_date(out_dir: Path = DEFAULT_OUT_DIR) -> str:
    """Try to read the simulation start date from the configured params file."""
    params_paths = []

    env_params = os.environ.get("CGEM_PARAMS_PATH")
    if env_params:
        params_paths.append(Path(env_params))

    params_paths.extend([
        Path("INPUT/params.txt"),
        Path("params.txt"),
        out_dir.parent / "INPUT/params.txt",
        out_dir.parent / "params.txt",
    ])

    seen = set()
    for p in params_paths:
        p = Path(p)
        if p in seen or not p.exists():
            continue
        seen.add(p)

        params = load_params(str(p))
        for key in ("SIM_START_DATE", "sim_start_date", "start_date", "START_DATE"):
            value = params.get(key)
            if value:
                return value
    return "2017-01-01"  # fallback default


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_geometry(out_dir: Path = DEFAULT_OUT_DIR) -> pd.DataFrame:
    """Load the geometry snapshot from OUT/geometry.csv.

    Returns DataFrame with columns: Cell, Distance, RiverbedDepth, Width,
    BaseArea, Chezy, FRIC.  Distance is converted to km.
    """
    path = out_dir / "geometry.csv"
    if not path.exists():
        # Try INPUT/Geometry/Geometry.csv as fallback
        fallback = Path("INPUT/Geometry/Geometry.csv")
        if fallback.exists():
            df = pd.read_csv(fallback)
            df.columns = [c.strip() for c in df.columns]
            if "Location" in df.columns:
                df = df.rename(columns={"Location": "Distance_km"})
            else:
                df["Distance_km"] = df.index
            return df
        raise FileNotFoundError(
            f"No geometry file found. Run the model first to generate {path}, "
            f"or place geometry data in INPUT/Geometry/Geometry.csv"
        )

    df = pd.read_csv(path)
    df["Distance_km"] = df["Distance"] / 1000.0
    return df


def load_output(var_name: str, out_dir: Path = DEFAULT_OUT_DIR,
                start_date: Optional[str] = None) -> pd.DataFrame:
    """Load a model output variable as a DataFrame.

    Parameters
    ----------
    var_name : str
        Variable name (e.g. 'Sal', 'O2', 'Reaction_NPP_NO3')
    out_dir : Path
        Output directory (default: OUT/)
    start_date : str, optional
        Simulation start date for binary files (default: read from params.txt)

    Returns
    -------
    pd.DataFrame
        Columns are cell distances (km), index is datetime.
    """
    csv_path = out_dir / f"{var_name}.csv"
    bin_path = out_dir / f"{var_name}.bin"

    if csv_path.exists():
        df = _read_csv(csv_path)
        # Convert cell column names to distance (km)
        geom = load_geometry(out_dir)
        if len(geom) >= len(df.columns):
            new_cols = geom["Distance_km"].values[:len(df.columns)]
            df.columns = new_cols
        return df

    elif bin_path.exists():
        cell_count, dt, dx, time_arr, values = _read_binary(bin_path)
        # Build datetime index
        if start_date is None:
            start_date = _detect_start_date(out_dir)
        start = pd.Timestamp(start_date)
        dates = [start + pd.Timedelta(seconds=float(t)) for t in time_arr]
        # Build distance columns
        distances_km = np.arange(cell_count) * dx / 1000.0
        df = pd.DataFrame(values, index=pd.DatetimeIndex(dates), columns=distances_km)
        df.index.name = "Date"
        return df

    else:
        available = list_variables(out_dir)
        raise FileNotFoundError(
            f"Variable '{var_name}' not found in {out_dir}.\n"
            f"Available variables: {', '.join(sorted(available))}"
        )


def list_variables(out_dir: Path = DEFAULT_OUT_DIR) -> list:
    """List all available output variables in the output directory."""
    out_dir = Path(out_dir)
    if not out_dir.exists():
        return []

    variables = set()
    for f in out_dir.iterdir():
        if f.suffix in (".csv", ".bin") and f.stem != "geometry":
            variables.add(f.stem)
    return sorted(variables)


def get_variable_info(var_name: str) -> tuple:
    """Return (long_name, unit) for a variable, or defaults if unknown."""
    return VARIABLE_INFO.get(var_name, (var_name, ""))


def load_params(params_path: str = "INPUT/params.txt") -> dict:
    """Parse params.txt into a dictionary of key-value pairs."""
    params = {}
    path = Path(params_path)
    if not path.exists():
        return params

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, _, val = line.partition("=")
                key = key.strip()
                val = val.split("#")[0].strip().strip('"')
                params[key] = val
    return params
