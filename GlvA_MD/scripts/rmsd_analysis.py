#!/usr/bin/env python3
"""
Compute RMSD for every group in a GROMACS index file.
Generates overlay plot and individual subplot grid.
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import parse_ndx, _to_ns

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
except ImportError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="RMSD analysis for all index groups")
    parser.add_argument("--tpr", required=True, help="GROMACS topology (.tpr)")
    parser.add_argument("--traj", required=True, help="Trajectory file (.xtc)")
    parser.add_argument("--ndx", required=True, help="Index file (.ndx)")
    parser.add_argument("--skip", nargs="*", default=[],
                        help="Group names to skip (space-separated)")
    parser.add_argument("--outdir", default="results", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    skip_set = set(args.skip) if args.skip else set()

    ndx_groups = parse_ndx(args.ndx)
    print(f"Loaded {len(ndx_groups)} groups from {args.ndx}")

    u = mda.Universe(args.tpr, args.traj)
    print(f"Loaded Universe with {len(u.trajectory)} frames")

    rmsd_all = {}
    for name, indices in ndx_groups.items():
        if name in skip_set:
            print(f"Skipping '{name}' (in skip list)...")
            continue
        sel = "index " + " ".join(map(str, indices))
        print(f"Calculating RMSD for {name} ({len(indices)} atoms)...")
        R = rms.RMSD(u, u, select=sel, superposition_selection="protein", ref_frame=0)
        R.run()
        rmsd_all[name] = R.results.rmsd

    if not rmsd_all:
        print("No groups computed.", file=sys.stderr)
        return

    print(f"\nDone. Computed RMSD for {len(rmsd_all)} groups.")

    data_rows = ["name,time_ps,rmsd"]
    for name, data in rmsd_all.items():
        for row in data:
            data_rows.append(f"{name},{row[1]},{row[2]}")
    csv_path = os.path.join(args.outdir, "rmsd_data.csv")
    with open(csv_path, "w") as f:
        f.write("\n".join(data_rows))
    print(f"Saved RMSD data to {csv_path}")

    plt.figure(figsize=(12, 7))
    for name, data in rmsd_all.items():
        plt.plot(_to_ns(data), data[:, 2], label=name, linewidth=1.2)
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (A)")
    plt.title("RMSD vs time for every index group")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=9)
    plt.tight_layout()
    overlay_path = os.path.join(args.outdir, "rmsd_overlay.png")
    plt.savefig(overlay_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved overlay plot to {overlay_path}")

    names = list(rmsd_all.keys())
    n = len(names)
    cols = 2
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(14, 3.2 * rows), squeeze=False)
    for i, name in enumerate(names):
        ax = axes[i // cols][i % cols]
        ax.plot(_to_ns(rmsd_all[name]), rmsd_all[name][:, 2], color="navy", linewidth=1.2)
        ax.set_title(name, fontsize=10)
        ax.set_ylabel("RMSD (A)")
        ax.set_xlabel("Time (ns)")
        ax.grid(True, linestyle="--", alpha=0.5)
    for j in range(n, rows * cols):
        axes[j // cols][j % cols].axis("off")
    plt.tight_layout()
    subplot_path = os.path.join(args.outdir, "rmsd_subplots.png")
    plt.savefig(subplot_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved subplot grid to {subplot_path}")


if __name__ == "__main__":
    main()
