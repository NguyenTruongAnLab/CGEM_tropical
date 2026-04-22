#!/usr/bin/env python3
"""
check_inputs.py — Comprehensive input data quality validation with diagnostic figure.

Scans all CSV boundary-condition and tributary files in INPUT/ for:
  1. File existence and readability
  2. Required columns (Time, variable columns)
  3. Temporal continuity (gaps, duplicates)
  4. Physical range plausibility
  5. Coefficient of variation (flags static/noisy signals)
  6. NaN percentage

Produces ONE comprehensive diagnostic figure with time series panels
for all input files, annotated with gap markers and range bands.

Usage:
    python tools/check_inputs.py
    python tools/check_inputs.py --input-dir INPUT --output-dir OUT

Outputs:
    OUT/figures/fig_input_diagnostics.png     (comprehensive figure)
    OUT/tables/check_inputs_report.csv        (machine-readable summary)
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cgem_read

# ---------------------------------------------------------------------------
# Plausible physical ranges (wide — flag gross errors, not calibration)
# ---------------------------------------------------------------------------

DEFAULT_RANGES = {
    "Discharge":     (0,    1e6,   "m³/s"),
    "Salinity":      (0,    40,    "PSU"),
    "Temperature":   (0,    45,    "°C"),
    "O2":            (0,    500,   "mmol/m³"),
    "NH4":           (0,    200,   "mmol/m³"),
    "NO3":           (0,    200,   "mmol/m³"),
    "PO4":           (0,    50,    "mmol/m³"),
    "SPM":           (0,    10,    "kg/m³"),
    "TOC":           (0,    2000,  "mmol/m³"),
    "DIC":           (0,    5000,  "mmol/m³"),
    "TAlk":          (0,    10000, "µeq/L"),
    "AT":            (0,    10000, "µeq/L"),
    "Phy":           (0,    500,   "mmolC/m³"),
    "Si":            (0,    500,   "mmol/m³"),
    "WaterLevel":    (-10,  10,    "m"),
    "elevation":     (-10,  10,    "m"),
    "waterDepth":    (0,    50,    "m"),
    "Velocity":      (-5,   5,     "m/s"),
    "PAR":           (0,    1000,  "W/m²"),
    "Light":         (0,    1000,  "W/m²"),
    "Wind":          (0,    30,    "m/s"),
    "pCO2":          (0,    50000, "µatm"),
}


# ---------------------------------------------------------------------------
# Quality check engine
# ---------------------------------------------------------------------------

def _iter_input_csv_files(input_dir):
    """Yield forcing-related CSV files that matter for a newcomer run."""
    candidate_dirs = [
        input_dir / "Boundary",
        input_dir / "Tributaries",
        input_dir / "carbonate",
    ]

    csv_files = []
    for folder in candidate_dirs:
        if folder.exists():
            csv_files.extend(
                sorted(
                    f for f in folder.rglob("*.csv")
                    if "summary" not in f.stem.lower()
                )
            )
    return csv_files

def _looks_like_headerless_two_column(filepath):
    """Detect numeric time,value files without a header row."""
    splitter = re.compile(r"[\s,]+")
    try:
        with open(filepath, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#") or line.startswith("%"):
                    continue
                parts = [token for token in splitter.split(line) if token]
                if len(parts) != 2:
                    return False
                try:
                    float(parts[0])
                    float(parts[1])
                    return True
                except ValueError:
                    return False
    except OSError:
        return False
    return False


def _read_input_csv(filepath):
    """Read either headerless numeric forcing files or regular CSV tables."""
    if _looks_like_headerless_two_column(filepath):
        return pd.read_csv(
            filepath,
            comment="#",
            header=None,
            names=["Time", "Value"],
            sep=r"[\s,]+",
            engine="python",
        )

    return pd.read_csv(filepath, comment="#")

def _find_time_col(df):
    """Find the time/date column in a DataFrame."""
    for candidate in ("Time", "time", "DATE", "Date", "datetime"):
        if candidate in df.columns:
            return candidate
    return None


def _get_range(col_name, context=""):
    """Match a column name to its plausible range."""
    labels = [str(col_name)]
    if context:
        labels.append(str(context))

    for keyword, (lo, hi, unit) in DEFAULT_RANGES.items():
        lowered = keyword.lower()
        for label in labels:
            tokens = [token for token in re.split(r"[^a-zA-Z0-9]+", label.lower()) if token]
            if label.lower() == lowered:
                return lo, hi, unit
            if lowered in tokens:
                return lo, hi, unit
            if len(lowered) >= 3 and any(token.startswith(lowered) for token in tokens):
                return lo, hi, unit
    return None, None, None


def _check_csv(filepath, report):
    """Check a single CSV file for issues."""
    name = filepath.name
    try:
        df = _read_input_csv(filepath)
    except Exception as e:
        report.append({"File": name, "Check": "readability", "Status": "FAIL",
                        "Detail": str(e)})
        return None

    report.append({"File": name, "Check": "readability", "Status": "PASS",
                    "Detail": f"{len(df)} rows × {len(df.columns)} columns"})

    time_col = _find_time_col(df)
    if time_col:
        parsed_time = pd.to_datetime(df[time_col], errors="coerce")
        if parsed_time.notna().sum() >= 2:
            times = pd.Series(parsed_time)
            n_dup = times.duplicated().sum()
            if n_dup > 0:
                report.append({"File": name, "Check": "duplicates", "Status": "WARN",
                                "Detail": f"{n_dup} duplicate timestamps"})
            else:
                report.append({"File": name, "Check": "duplicates", "Status": "PASS",
                                "Detail": "no duplicates"})

            dt = times.diff().dropna()
            if len(dt) > 1:
                median_dt = dt.median()
                gaps = dt[dt > 2 * median_dt]
                if len(gaps) > 0:
                    report.append({"File": name, "Check": "continuity", "Status": "WARN",
                                    "Detail": f"{len(gaps)} gaps > 2× median timestep"})
                else:
                    report.append({"File": name, "Check": "continuity", "Status": "PASS",
                                    "Detail": f"median dt = {median_dt}"})
        else:
            numeric_time = pd.to_numeric(df[time_col], errors="coerce")
            if numeric_time.notna().sum() >= 2:
                times = pd.Series(numeric_time)
                n_dup = times.duplicated().sum()
                if n_dup > 0:
                    report.append({"File": name, "Check": "duplicates", "Status": "WARN",
                                    "Detail": f"{n_dup} duplicate time indices"})
                else:
                    report.append({"File": name, "Check": "duplicates", "Status": "PASS",
                                    "Detail": "no duplicates"})

                dt = times.diff().dropna()
                if len(dt) > 1:
                    median_dt = dt.median()
                    gaps = dt[dt > 2 * median_dt]
                    if len(gaps) > 0:
                        report.append({"File": name, "Check": "continuity", "Status": "WARN",
                                        "Detail": f"{len(gaps)} gaps > 2× median timestep"})
                    else:
                        report.append({"File": name, "Check": "continuity", "Status": "PASS",
                                        "Detail": f"median step = {median_dt:g}"})
            else:
                report.append({"File": name, "Check": "time_parse", "Status": "WARN",
                                "Detail": f"Could not parse '{time_col}' as datetime or numeric time"})

    numeric_cols = [c for c in df.select_dtypes(include=[np.number]).columns if c != time_col]
    for col in numeric_cols:
        vals = df[col].dropna()
        if len(vals) == 0:
            report.append({"File": name, "Check": f"values:{col}", "Status": "WARN",
                            "Detail": "all NaN"})
            continue

        n_nan = df[col].isna().sum()
        if n_nan > 0:
            pct = 100.0 * n_nan / len(df)
            status = "WARN" if pct > 5 else "PASS"
            report.append({"File": name, "Check": f"missing:{col}", "Status": status,
                            "Detail": f"{n_nan} NaN ({pct:.1f}%)"})

        vmin, vmax = vals.min(), vals.max()
        lo, hi, unit = _get_range(col, filepath.stem)
        if lo is not None:
            if vmin < lo or vmax > hi:
                report.append({"File": name, "Check": f"range:{col}", "Status": "WARN",
                                "Detail": f"[{vmin:.3g}, {vmax:.3g}] outside [{lo}, {hi}] {unit}"})
            else:
                report.append({"File": name, "Check": f"range:{col}", "Status": "PASS",
                                "Detail": f"[{vmin:.3g}, {vmax:.3g}] within [{lo}, {hi}]"})

        if abs(vals.mean()) > 1e-15:
            cv = vals.std() / abs(vals.mean())
            if cv < 0.001:
                report.append({"File": name, "Check": f"variability:{col}", "Status": "WARN",
                                "Detail": f"CV = {cv:.6f} (nearly constant)"})

    return df


# ---------------------------------------------------------------------------
# Diagnostic figure
# ---------------------------------------------------------------------------

def _plot_diagnostic_figure(input_dir, fig_path, report_df):
    """Create comprehensive input-data diagnostic figure."""
    csv_files = _iter_input_csv_files(input_dir)
    if not csv_files:
        return

    # Collect plottable panels: (file_name, time_series, col_name, range_info)
    panels = []
    for fpath in csv_files:
        try:
            df = _read_input_csv(fpath)
        except Exception:
            continue
        time_col = _find_time_col(df)
        if not time_col:
            continue
        parsed_datetime = pd.to_datetime(df[time_col], errors="coerce")
        if parsed_datetime.notna().sum() >= 2:
            times = pd.Series(parsed_datetime)
            is_datetime = True
        else:
            numeric_time = pd.to_numeric(df[time_col], errors="coerce")
            if numeric_time.notna().sum() < 2:
                continue
            times = pd.Series(numeric_time)
            is_datetime = False

        numeric_cols = [c for c in df.select_dtypes(include=[np.number]).columns if c != time_col]
        # Limit to first 4 numeric cols per file to avoid plot overload
        for col in list(numeric_cols)[:4]:
            lo, hi, unit = _get_range(col, fpath.stem)
            panels.append({
                "file": fpath.stem,
                "col": col,
                "times": times,
                "values": pd.to_numeric(df[col], errors="coerce"),
                "lo": lo, "hi": hi, "unit": unit or "",
                "is_datetime": is_datetime,
            })

    if not panels:
        return

    # Layout: grid with max 3 columns
    n = len(panels)
    n_cols = min(3, n)
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(6.5 * n_cols, 3.0 * n_rows),
                             squeeze=False)

    for idx, p in enumerate(panels):
        ax = axes[idx // n_cols][idx % n_cols]
        times = p["times"]
        vals = p["values"]
        valid = vals.notna()

        # Plot time series
        ax.plot(times[valid], vals[valid], "b-", lw=0.6, alpha=0.8)

        # Range band
        if p["lo"] is not None:
            ax.axhspan(p["lo"], p["hi"], color="green", alpha=0.06)
            # Mark out-of-range points
            oor = valid & ((vals < p["lo"]) | (vals > p["hi"]))
            if oor.any():
                ax.scatter(times[oor], vals[oor], c="red", s=8, zorder=5,
                           label=f"{oor.sum()} out-of-range")
                ax.legend(fontsize=7, loc="upper right")

        # Mark gaps
        dt = times.diff()
        if len(dt) > 1:
            median_dt = dt.dropna().median()
            gap_idx = dt > 2 * median_dt
            if gap_idx.any():
                for t in times[gap_idx]:
                    ax.axvline(t, color="orange", lw=0.5, alpha=0.5)

        # NaN fraction annotation
        nan_pct = 100.0 * (~valid).sum() / len(vals)
        status_color = "red" if nan_pct > 10 else ("orange" if nan_pct > 1 else "green")
        ax.text(0.98, 0.95, f"NaN: {nan_pct:.1f}%",
                transform=ax.transAxes, ha="right", va="top", fontsize=8,
                color=status_color, fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7))

        # Title and labels
        title = f"{p['file']}: {p['col']}"
        if p["unit"]:
            title += f" ({p['unit']})"
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.grid(True, ls=":", alpha=0.3)
        if p["is_datetime"]:
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
            ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=3, maxticks=6))
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, fontsize=7, ha="right")
            ax.set_xlabel("Time")
        else:
            ax.set_xlabel("Time index")
            plt.setp(ax.xaxis.get_majorticklabels(), fontsize=7)

    # Hide unused
    for idx in range(n, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)

    plt.tight_layout(pad=2.0)
    fig.savefig(fig_path, dpi=250, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def check(input_dir=None, out_dir=None):
    """Scan all CSV files in input_dir and produce quality report + figure."""
    input_dir = Path(input_dir or "INPUT")
    out_dir = Path(out_dir or "OUT")
    fig_dir = out_dir / "figures"
    tbl_dir = out_dir / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir.mkdir(parents=True, exist_ok=True)

    csv_files = _iter_input_csv_files(input_dir)
    if not csv_files:
        print(f"  No CSV files found in {input_dir}")
        return pd.DataFrame()

    report = []
    for f in csv_files:
        _check_csv(f, report)

    df = pd.DataFrame(report)
    report_path = tbl_dir / "check_inputs_report.csv"
    df.to_csv(report_path, index=False)

    # Diagnostic figure
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fig_path = fig_dir / "fig_input_diagnostics.png"
    _plot_diagnostic_figure(input_dir, fig_path, df)

    # Console summary
    n_pass = (df["Status"] == "PASS").sum()
    n_warn = (df["Status"] == "WARN").sum()
    n_fail = (df["Status"] == "FAIL").sum()

    print(f"\n{'='*65}")
    print("  INPUT DATA QUALITY REPORT")
    print(f"{'='*65}")
    print(f"  Files scanned : {len(csv_files)}")
    print(f"  Checks        : {len(df)}  (PASS={n_pass}  WARN={n_warn}  FAIL={n_fail})")
    print(f"{'='*65}")

    if n_fail > 0:
        print("\n  FAILURES:")
        for _, row in df[df["Status"] == "FAIL"].iterrows():
            print(f"    X {row['File']}: {row['Check']} -- {row['Detail']}")

    if n_warn > 0:
        print(f"\n  WARNINGS ({n_warn}):")
        for _, row in df[df["Status"] == "WARN"].iterrows():
            print(f"    ! {row['File']}: {row['Check']} -- {row['Detail']}")

    print(f"\n  Report  -> {report_path}")
    if fig_path.exists():
        fig_path_ts = fig_dir / f"fig_input_diagnostics_{ts}.png"
        import shutil
        shutil.copy2(fig_path, fig_path_ts)
        print(f"  Figure  -> {fig_path}")
    print(f"{'='*65}")
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate model input CSV files")
    parser.add_argument("--input-dir", default="INPUT", help="Input directory to scan")
    parser.add_argument("--output-dir", default="OUT", help="Where to save report")
    args = parser.parse_args()
    check(args.input_dir, args.output_dir)
