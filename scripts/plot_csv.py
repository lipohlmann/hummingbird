#!/usr/bin/env python3
"""Plot flux solutions from a hummingbird simulation results CSV (see
Results::ToCSV in src/utils/results.cc for the file format: node_id,x,y,z,
scalar_flux, then angular_flux_i,azimuth_i,polar_i per ordinate i)."""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

DPI = 300
FIGSIZE = (6.0, 4.0)
POSITION_TOLERANCE = 1e-9
X_LABEL = "x [cm]"
SCALAR_FLUX_LABEL = "Scalar Flux [n/cm²/s]"
ANGULAR_FLUX_LABEL = "Angular Flux [n/cm²/s/str]"


def configure_style():
    """Serifed, publication-style figure defaults."""
    sns.set_theme(style="ticks")
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Nimbus Roman", "Liberation Serif", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "legend.fontsize": 10.5,
        "legend.title_fontsize": 11,
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "legend.edgecolor": "0.8",
        "axes.linewidth": 1.0,
        "lines.linewidth": 1.8,
        "lines.markersize": 5,
        "lines.markeredgewidth": 0,
    })


def finalize_axes(fig, ax):
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.4)
    sns.despine(fig=fig)
    fig.tight_layout()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot flux solutions from a hummingbird simulation results CSV.")
    parser.add_argument("csv_file", type=Path,
                        help="Path to a *_results.csv file produced by the simulation.")
    parser.add_argument("--angular", action="store_true",
                        help="Also produce a single combined plot of every angular flux "
                             "value, with a legend differentiating each ordinate's angle. "
                             "By default, only the scalar flux is plotted.")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Directory to save plots in (default: the CSV file's directory).")
    return parser.parse_args()


def check_is_1d(df):
    for axis in ("y", "z"):
        if (df[axis] - df[axis].iloc[0]).abs().max() > POSITION_TOLERANCE:
            raise ValueError(
                f"CSV data is not 1D: the '{axis}' coordinate varies across nodes. "
                "This script only supports plotting 1D simulation output.")


def count_ordinates(df):
    n = 0
    while f"angular_flux_{n}" in df.columns:
        n += 1
    return n


def plot_scalar_flux(df, palette, output_path):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    sns.lineplot(data=df, x="x", y="scalar_flux", marker="o", ax=ax, color=palette[0])
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(SCALAR_FLUX_LABEL)
    ax.set_title("Scalar Flux")
    finalize_axes(fig, ax)
    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_angular_fluxes(df, n_ordinates, palette, output_path):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for i in range(n_ordinates):
        azimuth = df[f"azimuth_{i}"].iloc[0]
        polar = df[f"polar_{i}"].iloc[0]
        label = f"azimuth={azimuth:.3f} rad, polar={polar:.3f} rad"
        sns.lineplot(data=df, x="x", y=f"angular_flux_{i}", marker="o", ax=ax,
                    color=palette[i % len(palette)], label=label)
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(ANGULAR_FLUX_LABEL)
    ax.set_title("Angular Flux by Ordinate")
    ax.legend(title="Ordinate")
    finalize_axes(fig, ax)
    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()

    df = pd.read_csv(args.csv_file).sort_values("x")
    check_is_1d(df)

    output_dir = args.output_dir or args.csv_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.csv_file.stem

    configure_style()
    palette = sns.color_palette("colorblind")

    plot_scalar_flux(df, palette, output_dir / f"{stem}_scalar_flux.png")

    if args.angular:
        n_ordinates = count_ordinates(df)
        if n_ordinates == 0:
            raise ValueError(
                "--angular was given but the CSV has no angular_flux_* columns.")
        plot_angular_fluxes(df, n_ordinates, palette,
                            output_dir / f"{stem}_angular_flux.png")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
