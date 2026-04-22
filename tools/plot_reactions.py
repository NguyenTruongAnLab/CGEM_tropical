#!/usr/bin/env python3
"""
plot_reactions.py — Comprehensive biogeochemical reaction-rate visualization.

Produces ONE multi-panel figure showing longitudinal profiles of all
reaction rates grouped by process: production, nutrient cycling,
gas exchange, carbonate system, limitation factors, and light attenuation.
Each group shows temporal mean ± std with dry/wet seasonal overlays.

Requires enable_reaction_output = 1 in params.txt.

Usage:
    python tools/plot_reactions.py                   # All groups
    python tools/plot_reactions.py --mode limitations # Single group
    python tools/plot_reactions.py --snapshot         # Last timestep

Output:
    OUT/figures/fig_reactions.png
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cgem_read

RATE_GROUPS = {
    "production": {
        "title": "Primary Production & Mortality",
        "vars": ["Reaction_NPP_NO3", "Reaction_NPP_NH4", "Reaction_phydeath"],
        "colors": ["forestgreen", "limegreen", "indianred"],
    },
    "nutrients": {
        "title": "Nutrient Cycling",
        "vars": ["Reaction_nitrification", "Reaction_denitrification", "Reaction_adegradation"],
        "colors": ["darkorange", "purple", "saddlebrown"],
    },
    "gas_exchange": {
        "title": "Air-Water Gas Exchange",
        "vars": ["Reaction_O2_exchange", "Reaction_CO2_exchange"],
        "colors": ["steelblue", "gray"],
    },
    "carbonate": {
        "title": "Carbonate System Rates",
        "vars": ["Reaction_DIC", "Reaction_TA"],
        "colors": ["teal", "olive"],
    },
    "limitations": {
        "title": "Growth Limitation Factors",
        "vars": ["Diag_fN", "Diag_fP", "Diag_fSi", "Diag_fI"],
        "colors": ["darkorange", "purple", "goldenrod", "steelblue"],
    },
    "light": {
        "title": "Light Attenuation",
        "vars": ["Diag_KD"],
        "colors": ["saddlebrown"],
    },
}


def plot(out_dir, mode="all", snapshot=False, save_dir=None):
    """Create comprehensive reaction-rate diagnostic figure."""
    out_dir = Path(out_dir)
    if save_dir is None:
        save_dir = out_dir
    fig_dir = Path(save_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    available = cgem_read.list_variables(out_dir)

    if mode == "all":
        groups = [(k, v) for k, v in RATE_GROUPS.items()
                  if any(var in available for var in v["vars"])]
    else:
        if mode not in RATE_GROUPS:
            print(f"  Unknown mode '{mode}'")
            return
        g = RATE_GROUPS[mode]
        if not any(v in available for v in g["vars"]):
            print(f"  No data for group '{mode}'. Set enable_reaction_output=1.")
            return
        groups = [(mode, g)]

    if not groups:
        print("  No reaction rate or diagnostic output found.\n"
              "  Set enable_reaction_output = 1 in params.txt and re-run.")
        return

    n = len(groups)
    fig, axes = plt.subplots(n, 1, figsize=(12, 3.5 * n), sharex=True, squeeze=False)

    for row_idx, (gname, ginfo) in enumerate(groups):
        ax = axes[row_idx, 0]
        vars_to_plot = [v for v in ginfo["vars"] if v in available]
        last_unit = ""

        for var, color in zip(vars_to_plot, ginfo["colors"]):
            df = cgem_read.load_output(var, out_dir)
            distances = df.columns.astype(float).values
            long_name, unit = cgem_read.get_variable_info(var)
            if unit:
                last_unit = unit

            if snapshot:
                ax.plot(distances, df.iloc[-1].values, color=color,
                        lw=1.3, label=long_name)
            else:
                # Overall mean ± std
                mean_p = df.mean(axis=0).values
                std_p = df.std(axis=0).values
                ax.plot(distances, mean_p, color=color, lw=1.8, label=long_name)
                ax.fill_between(distances, mean_p - std_p, mean_p + std_p,
                                color=color, alpha=0.12)

                # Seasonal overlays (dashed)
                seasons = cgem_read.classify_season(df.index)
                for sn, ls in [("dry", ":"), ("wet", "--")]:
                    mask = seasons == sn
                    if mask.any():
                        ax.plot(distances, df.loc[mask].mean(axis=0).values,
                                color=color, lw=0.9, ls=ls, alpha=0.6)

        ax.set_title(ginfo["title"], fontsize=10, fontweight="bold")
        ax.axhline(0, color="gray", lw=0.5, ls="--")
        ax.grid(True, alpha=0.3)
        if last_unit:
            ax.set_ylabel(f"Rate ({last_unit})", fontsize=9)
        ax.legend(fontsize=8, loc="best", frameon=True)

    axes[-1, 0].set_xlabel("Distance from mouth (km)", fontsize=10)

    plt.tight_layout(pad=2.0)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fig_path = fig_dir / "fig_reactions.png"
    fig.savefig(fig_path, dpi=250, bbox_inches="tight")
    fig.savefig(fig_dir / f"fig_reactions_{ts}.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    print(f"\n{'='*65}")
    print("  BIOGEOCHEMICAL REACTION RATES")
    print(f"{'='*65}")
    for gname, ginfo in groups:
        vars_found = [v for v in ginfo["vars"] if v in available]
        print(f"  {ginfo['title']:35s}  [{', '.join(vars_found)}]")
    print(f"  Figure -> {fig_path}")
    print(f"{'='*65}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Biogeochemical reaction rates")
    parser.add_argument("--output-dir", default="OUT", help="Model output directory")
    parser.add_argument("--mode", default="all",
                        choices=["all", "production", "nutrients", "gas_exchange",
                                 "carbonate", "limitations", "light"],
                        help="Which rate group to plot")
    parser.add_argument("--snapshot", action="store_true",
                        help="Plot last timestep instead of time-mean")
    parser.add_argument("--save-dir", default=None, help="Output directory for figures")
    args = parser.parse_args()
    plot(Path(args.output_dir), args.mode, args.snapshot,
         args.save_dir or args.output_dir)
