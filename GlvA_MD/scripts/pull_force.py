#!/usr/bin/env python3
"""
Analysis of COM pulling force from GROMACS pullf.xvg output.

Reads the pullf XVG file, computes force statistics, applies a
Savitzky-Golay filter to smooth the noisy force curve, optionally
estimates mechanical work, and generates a plot.

Usage:
    python3 scripts/pull_force.py \
        --pullf step5_1_pullf.xvg \
        --smooth-window 51 \
        --smooth-polyorder 3 \
        --pulling-velocity 0.01 \
        --outdir results
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.signal import savgol_filter
except ImportError:
    print("ERROR: scipy is required (pip install scipy)", file=sys.stderr)
    sys.exit(1)


# ============================================================================
# XVG parser
# ============================================================================

def parse_xvg_force(filename):
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    time = []
    force = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(("@", "#")) or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                time.append(float(parts[0]))
                force.append(float(parts[1]))

    return np.array(time), np.array(force)


# ============================================================================
# Analysis
# ============================================================================

def analyze_force(time_ps, force_kj_mol_nm, smooth_window=51,
                  smooth_polyorder=3, pulling_velocity_nm_ps=None):
    mean_force = np.mean(force_kj_mol_nm)
    max_force = np.max(force_kj_mol_nm)
    min_force = np.min(force_kj_mol_nm)
    std_force = np.std(force_kj_mol_nm)
    time_at_max = time_ps[np.argmax(force_kj_mol_nm)]

    # Savitzky-Golay smoothing
    win = min(smooth_window,
              len(force_kj_mol_nm) - 1 if len(force_kj_mol_nm) % 2 != 0
              else len(force_kj_mol_nm) - 2)
    if win % 2 == 0:
        win -= 1
    if win < smooth_polyorder + 2:
        win = smooth_polyorder + 2
    if win % 2 == 0:
        win += 1

    smoothed = savgol_filter(force_kj_mol_nm, window_length=win,
                             polyorder=smooth_polyorder)

    work_kj_mol = None
    if pulling_velocity_nm_ps is not None:
        dt = np.diff(time_ps)
        f_mid = (force_kj_mol_nm[:-1] + force_kj_mol_nm[1:]) / 2
        dx = pulling_velocity_nm_ps * dt
        work_kj_mol = np.sum(f_mid * dx)

    return {
        "time_ps": time_ps,
        "force_kj_mol_nm": force_kj_mol_nm,
        "smoothed_force_kj_mol_nm": smoothed,
        "mean_force_kj_mol_nm": mean_force,
        "max_force_kj_mol_nm": max_force,
        "min_force_kj_mol_nm": min_force,
        "std_force_kj_mol_nm": std_force,
        "time_at_max_ps": time_at_max,
        "work_kj_mol": work_kj_mol,
    }


# ============================================================================
# Plot
# ============================================================================

def plot_results(results, outdir="results"):
    os.makedirs(outdir, exist_ok=True)

    plt.figure(figsize=(10, 6))
    plt.plot(results["time_ps"], results["force_kj_mol_nm"],
             color="lightgray", alpha=0.7, label="Raw force (noisy)")
    plt.plot(results["time_ps"], results["smoothed_force_kj_mol_nm"],
             color="blue", linewidth=2, label="Smoothed (Savitzky-Golay)")
    plt.axhline(results["max_force_kj_mol_nm"], color="red",
                linestyle="--", linewidth=1,
                label=f"Max force: {results['max_force_kj_mol_nm']:.2f} kJ/(mol·nm)")
    plt.axvline(results["time_at_max_ps"], color="red",
                linestyle=":", linewidth=1)

    plt.title("COM Pulling Force Analysis", fontsize=14)
    plt.xlabel("Time (ps)", fontsize=12)
    plt.ylabel("Force (kJ/(mol·nm))", fontsize=12)
    plt.legend(loc="best")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    outpath = os.path.join(outdir, "pull_force.png")
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot saved to {outpath}")


# ============================================================================
# Report
# ============================================================================

def print_report(results, pulling_velocity_nm_ps):
    sep = "=" * 60
    print(f"\n{sep}")
    print("COM PULLING FORCE ANALYSIS")
    print(sep)
    print(f"  Data points:             {len(results['time_ps'])}")
    print(f"  Duration:                "
          f"{results['time_ps'][-1] - results['time_ps'][0]:.2f} ps")
    print(f"  Mean force:              "
          f"{results['mean_force_kj_mol_nm']:.2f} kJ/(mol·nm)")
    print(f"  Max force:               "
          f"{results['max_force_kj_mol_nm']:.2f} kJ/(mol·nm)  "
          f"(at t = {results['time_at_max_ps']:.2f} ps)")
    print(f"  Min force:               "
          f"{results['min_force_kj_mol_nm']:.2f} kJ/(mol·nm)")
    print(f"  Std deviation (noise):   "
          f"{results['std_force_kj_mol_nm']:.2f} kJ/(mol·nm)")

    if results["work_kj_mol"] is not None:
        print(f"  Estimated work:          "
              f"{results['work_kj_mol']:.4f} kJ/mol  "
              f"(v = {pulling_velocity_nm_ps} nm/ps)")
    else:
        print("  Work estimation:         skipped "
              "(no pulling velocity provided)")

    print(sep)


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="COM pulling force analysis from pullf.xvg"
    )
    parser.add_argument("--pullf", default="step5_1_pullf.xvg",
                        help="Path to pullf.xvg file")
    parser.add_argument("--smooth-window", type=int, default=51,
                        help="Savitzky-Golay window size (odd)")
    parser.add_argument("--smooth-polyorder", type=int, default=3,
                        help="Savitzky-Golay polynomial order")
    parser.add_argument("--pulling-velocity", type=float, default=None,
                        help="Pulling velocity in nm/ps (for work estimation)")
    parser.add_argument("--outdir", default="results",
                        help="Output directory")
    args = parser.parse_args()

    if not os.path.isfile(args.pullf):
        print(f"ERROR: pullf file not found: {args.pullf}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading: {args.pullf}")
    time_ps, force = parse_xvg_force(args.pullf)

    results = analyze_force(
        time_ps, force,
        smooth_window=args.smooth_window,
        smooth_polyorder=args.smooth_polyorder,
        pulling_velocity_nm_ps=args.pulling_velocity,
    )
    print_report(results, args.pulling_velocity)
    plot_results(results, outdir=args.outdir)


if __name__ == "__main__":
    main()
