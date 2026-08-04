#!/usr/bin/env python3
"""
Calculate distances over time for specified atom pairs.
Generates a combined distance plot and CSV data.
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tqdm import tqdm

try:
    import MDAnalysis as mda
except ImportError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Bond distance analysis")
    parser.add_argument("--tpr", required=True, help="GROMACS topology (.tpr)")
    parser.add_argument("--traj", required=True, help="Trajectory file (.xtc)")
    parser.add_argument("--bonds", nargs="+", action="append",
                        help="Atom pair indices: --bonds IDX1 IDX2 [--bonds IDX3 IDX4 ...]")
    parser.add_argument("--labels", nargs="*", default=[],
                        help="Labels for each bond pair (same order as --bonds)")
    parser.add_argument("--outdir", default="results", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if not args.bonds:
        print("ERROR: At least one --bonds pair required", file=sys.stderr)
        sys.exit(1)

    bond_pairs = [(int(p[0]), int(p[1])) for p in args.bonds]

    labels = args.labels if args.labels else [f"Bond_{i+1}" for i in range(len(bond_pairs))]
    while len(labels) < len(bond_pairs):
        labels.append(f"Bond_{len(labels)+1}")

    u = mda.Universe(args.tpr, args.traj)
    print(f"Loaded Universe with {len(u.trajectory)} frames")

    dists = {lbl: [] for lbl in labels}
    times = []

    print("Computing bond distances...")
    for ts in tqdm(u.trajectory, desc="Trajectory"):
        times.append(ts.time)
        for i, (a, b) in enumerate(bond_pairs):
            dists[labels[i]].append(
                np.linalg.norm(u.atoms[a].position - u.atoms[b].position)
            )

    times_ns = np.array(times) / 1000.0

    rows = ["time_ns," + ",".join(labels)]
    for t_idx in range(len(times_ns)):
        vals = [str(times_ns[t_idx])]
        for lbl in labels:
            vals.append(str(dists[lbl][t_idx]))
        rows.append(",".join(vals))
    csv_path = os.path.join(args.outdir, "bond_distances.csv")
    with open(csv_path, "w") as f:
        f.write("\n".join(rows))
    print(f"Saved bond distances to {csv_path}")

    fig, axes = plt.subplots(len(bond_pairs), 1, figsize=(12, 3 * len(bond_pairs)), sharex=True)
    if len(bond_pairs) == 1:
        axes = [axes]
    for i, lbl in enumerate(labels):
        axes[i].plot(times_ns, dists[lbl], color="firebrick", linewidth=1.2)
        axes[i].set_ylabel("Distance (A)")
        axes[i].set_title(lbl)
        axes[i].grid(True, linestyle="--", alpha=0.6)
    axes[-1].set_xlabel("Time (ns)")
    plt.tight_layout()
    plot_path = os.path.join(args.outdir, "bond_distances.png")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved bond distance plot to {plot_path}")


if __name__ == "__main__":
    main()
