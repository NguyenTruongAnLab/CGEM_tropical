#!/usr/bin/env python3
"""Convert field-observation CSVs to C-GEM model units.

By default this script reads the example observation tables in
`INPUT/Validation/` and writes converted copies to
`INPUT/Validation/converted/` so the source files are never overwritten.

The repository-shipped `obs_*.csv` files are already in model units; this
script is mainly intended for new field tables prepared in common lab units.
See `INPUT/Validation/README.md` for the expected target schema.
"""

import argparse
import shutil
from pathlib import Path

import pandas as pd

VAL_DIR = Path("INPUT/Validation")

# ── Field → Model conversion factors ──────────────────────────
# obs_value * factor = model_value
CONVERSIONS = {
    # profiles CSV columns → model name + factor
    "Salinity_PSU": ("Sal", 1.0),
    "DO_mgL": ("O2", 1000.0 / 32.00),
    "NH4_mgNL": ("NH4", 1000.0 / 14.01),
    "NO3_mgNL": ("NO3", 1000.0 / 14.01),
    "PO4_mgPL": ("PO4", 1000.0 / 30.97),
    "TOC_mgCL": ("TOC", 1000.0 / 12.01),
    "TSS_mgL": ("SPM", 0.001),  # mg/L → kg/m³
    "Chla_ugL": ("Phy1", 15.0 / 12.01),  # C:Chl=15
    "DSi_mgSiL": ("Si", 1000.0 / 28.09),
    "DIC_mgCL": ("DIC", 1000.0 / 12.01),
    "Alkalinity_mgCaCO3L": ("AT", 1000.0 / 50.04345),
    "pCO2_uatm": ("pCO2", 1.0),
    "pH": ("PH", 1.0),
}

# Columns copied as-is (metadata)
META_COLS = ["Date", "Station", "Distance_km", "Season"]


def convert_csv(src: Path, dst: Path, extra_meta=None):
    """Read field-unit CSV, convert to model units, and save a new file."""
    df = pd.read_csv(src)
    out = pd.DataFrame()

    # Copy metadata columns
    all_meta = META_COLS + (extra_meta or [])
    for col in all_meta:
        if col in df.columns:
            out[col] = df[col]

    # Convert data columns
    for field_col, (model_col, factor) in CONVERSIONS.items():
        if field_col in df.columns:
            out[model_col] = pd.to_numeric(df[field_col], errors="coerce") * factor

    out.to_csv(dst, index=False, float_format="%.6g")
    print(f"  {src.name} → {dst.name}  ({len(out)} rows, {len(out.columns)} cols)")
    return out


def _copy_if_present(src: Path, dst: Path, overwrite=False):
    if not src.exists():
        return False
    if dst.exists() and not overwrite:
        print(f"  Skipping {dst.name} (already exists; use --overwrite to replace)")
        return False
    shutil.copy2(src, dst)
    print(f"  {src.name} → {dst.name}  (copied; already in model units)")
    return True


def convert_default_tables(input_dir: Path, output_dir: Path, overwrite=False):
    """Convert the default example tables without modifying the originals."""
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Converting field-unit observation tables to C-GEM model units...")

    # 1. Longitudinal profiles
    profiles_src = input_dir / "obs_longitudinal_profiles.csv"
    profiles_dst = output_dir / "obs_longitudinal_profiles_model_units.csv"
    if profiles_src.exists():
        if profiles_dst.exists() and not overwrite:
            print(f"  Skipping {profiles_dst.name} (already exists; use --overwrite to replace)")
        else:
            convert_csv(profiles_src, profiles_dst)

    # 2. Time series
    ts_src = input_dir / "obs_timeseries.csv"
    ts_dst = output_dir / "obs_timeseries_model_units.csv"
    if ts_src.exists():
        if ts_dst.exists() and not overwrite:
            print(f"  Skipping {ts_dst.name} (already exists; use --overwrite to replace)")
        else:
            convert_csv(ts_src, ts_dst, extra_meta=["Source"])

    # 3. Tidal range (already in model/physical units, so copy only)
    tidal_src = input_dir / "obs_tidal_range.csv"
    tidal_dst = output_dir / "obs_tidal_range_model_units.csv"
    _copy_if_present(tidal_src, tidal_dst, overwrite=overwrite)

    print(f"Done. Converted files are in {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="Convert observation CSVs to C-GEM model units")
    parser.add_argument("--input-dir", default=str(VAL_DIR), help="Directory containing observation CSV files")
    parser.add_argument("--output-dir", default=str(VAL_DIR / "converted"),
                        help="Directory for converted CSV copies")
    parser.add_argument("--overwrite", action="store_true",
                        help="Replace converted files if they already exist")
    args = parser.parse_args()

    convert_default_tables(Path(args.input_dir), Path(args.output_dir), overwrite=args.overwrite)


if __name__ == "__main__":
    main()
