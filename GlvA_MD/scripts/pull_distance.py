#!/usr/bin/env python3
"""
Analysis of COM pulling distance from GROMACS pullx.xvg output.

Reads the pullx XVG file, computes thermal fluctuation statistics
(RMSF, effective stiffness), and generates a plot with mean distance,
±1σ band, and the flat-bottom potential boundary.

Usage:
    python3 scripts/pull_distance.py \
        --pullx step5_1_pullx.xvg \
        --temperature 300 \
        --flat-bottom 0.45 \
        --outdir results
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ============================================================================
# XVG parser
# ============================================================================

def parse_xvg_distance(filename):
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    time = []
    distance = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(("@", "#")) or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                time.append(float(parts[0]))
                distance.append(float(parts[1]))

    return np.array(time), np.array(distance)


# ============================================================================
# Analysis
# ============================================================================

def analyze_thermal_fluctuations(time_ps, distance_nm, temperature_k=300.0):
    kb = 0.008314462618  # kJ / (mol * K)
    kt = kb * temperature_k  # kJ / mol

    mean_dist = np.mean(distance_nm)
    std_dist = np.std(distance_nm)
    min_dist = np.min(distance_nm)
    max_dist = np.max(distance_nm)

    variance = std_dist ** 2
    k_eff = kt / variance if variance > 0 else 0.0

    return {
        "mean_dist_nm": mean_dist,
        "std_dist_nm": std_dist,
        "min_dist_nm": min_dist,
        "max_dist_nm": max_dist,
        "k_eff": k_eff,
        "kt": kt,
    }


# ============================================================================
# Plot
# ============================================================================

def plot_results(time_ps, distance_nm, results, flat_bottom_nm=0.45,
                 outdir="results"):
    os.makedirs(outdir, exist_ok=True)

    plt.figure(figsize=(10, 5))
    plt.plot(time_ps, distance_nm, color="blue", alpha=0.5,
             linewidth=0.5, label="Distance")
    plt.axhline(results["mean_dist_nm"], color="red", linestyle="--",
                label=f"Mean: {results['mean_dist_nm']:.3f} nm")
    plt.fill_between(time_ps,
                     results["mean_dist_nm"] - results["std_dist_nm"],
                     results["mean_dist_nm"] + results["std_dist_nm"],
                     color="red", alpha=0.2, label="±1σ (thermal fluctuations)")
    plt.axhline(flat_bottom_nm, color="green", linestyle=":", linewidth=2,
                label=f"Flat-bottom boundary ({flat_bottom_nm} nm)")

    plt.title("COM Pulling Distance — Thermal Fluctuations", fontsize=14)
    plt.xlabel("Time (ps)", fontsize=12)
    plt.ylabel("Distance (nm)", fontsize=12)
    plt.legend(loc="best")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    outpath = os.path.join(outdir, "pull_distance.png")
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot saved to {outpath}")


# ============================================================================
# Report
# ============================================================================

def print_report(results, flat_bottom_nm, temperature_k):
    sep = "=" * 60
    print(f"\n{sep}")
    print("COM PULLING DISTANCE — THERMAL FLUCTUATION ANALYSIS")
    print(sep)
    print(f"  Temperature:              {temperature_k} K  "
          f"(k_B T = {results['kt']:.4f} kJ/mol)")
    print(f"  Mean distance:            {results['mean_dist_nm']:.4f} nm")
    print(f"  Min / Max distance:       {results['min_dist_nm']:.4f}  /  "
          f"{results['max_dist_nm']:.4f} nm")
    print(f"  Fluctuation amplitude (RMSF): {results['std_dist_nm']:.4f} nm  "
          f"({results['std_dist_nm'] * 10:.2f} Å)")
    print(f"  Effective stiffness (k_eff):  {results['k_eff']:.1f} kJ/(mol·nm²)")
    print(f"  Flat-bottom boundary:     {flat_bottom_nm} nm")
    print(sep)


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="COM pulling distance analysis from pullx.xvg"
    )
    parser.add_argument("--pullx", default="step5_1_pullx.xvg",
                        help="Path to pullx.xvg file")
    parser.add_argument("--temperature", type=float, default=300.0,
                        help="Temperature in K")
    parser.add_argument("--flat-bottom", type=float, default=0.45,
                        help="Flat-bottom potential boundary (nm)")
    parser.add_argument("--outdir", default="results",
                        help="Output directory")
    args = parser.parse_args()

    if not os.path.isfile(args.pullx):
        print(f"ERROR: pullx file not found: {args.pullx}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading: {args.pullx}")
    time_ps, distance_nm = parse_xvg_distance(args.pullx)

    results = analyze_thermal_fluctuations(
        time_ps, distance_nm, temperature_k=args.temperature
    )
    print_report(results, args.flat_bottom, args.temperature)
    plot_results(time_ps, distance_nm, results,
                 flat_bottom_nm=args.flat_bottom, outdir=args.outdir)


if __name__ == "__main__":
    main()
