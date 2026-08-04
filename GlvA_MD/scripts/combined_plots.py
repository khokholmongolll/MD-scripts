#!/usr/bin/env python3
"""
Generate combined figure: RMSD plots + bond distance plots in a single vertical stack.
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Combined RMSD + bond distance plot")
    parser.add_argument("--outdir", default="results", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    rmsd_csv = os.path.join(args.outdir, "rmsd_data.csv")
    bond_csv = os.path.join(args.outdir, "bond_distances.csv")

    has_rmsd = os.path.isfile(rmsd_csv)
    has_bond = os.path.isfile(bond_csv)

    if not has_rmsd and not has_bond:
        print("No RMSD or bond distance data found. Run rmsd_analysis.py and "
              "bond_distance.py first.")
        return

    sections = []

    if has_rmsd:
        rmsd_data = {}
        with open(rmsd_csv, "r") as f:
            lines = f.readlines()
        header = lines[0].strip().split(",")
        for line in lines[1:]:
            name, time_ps, val = line.strip().split(",", 2)
            if name not in rmsd_data:
                rmsd_data[name] = {"time": [], "rmsd": []}
            rmsd_data[name]["time"].append(float(time_ps))
            rmsd_data[name]["rmsd"].append(float(val))
        sections.append(("rmsd", rmsd_data))

    if has_bond:
        with open(bond_csv, "r") as f:
            lines = f.readlines()
        header = lines[0].strip().split(",")
        bond_data = {h: {"time": [], "dist": []} for h in header[1:]}
        for line in lines[1:]:
            parts = line.strip().split(",")
            tn = float(parts[0])
            for i, lbl in enumerate(header[1:]):
                bond_data[lbl]["time"].append(tn)
                bond_data[lbl]["dist"].append(float(parts[i + 1]))
        sections.append(("bond", bond_data))

    n_plots = 0
    if has_rmsd:
        n_plots += len(sections[0][1])
    if has_bond:
        n_plots += len(sections[1][1])

    fig, axes = plt.subplots(n_plots, 1, figsize=(12, 3 * n_plots), sharex=False)
    if n_plots == 1:
        axes = [axes]

    idx = 0
    for stype, sec_data in sections:
        if stype == "rmsd":
            for name, d in sec_data.items():
                t_ns = np.array(d["time"]) / 1000.0
                axes[idx].plot(t_ns, d["rmsd"], color='navy')
                axes[idx].set_ylabel("RMSD (A)")
                axes[idx].set_title(name)
                axes[idx].grid(True, linestyle='--', alpha=0.6)
                idx += 1
        elif stype == "bond":
            for lbl, d in sec_data.items():
                axes[idx].plot(d["time"], d["dist"], color='firebrick')
                axes[idx].set_ylabel("Distance (A)")
                axes[idx].set_title(lbl)
                axes[idx].grid(True, linestyle='--', alpha=0.6)
                idx += 1

    axes[-1].set_xlabel("Time (ns)")
    plt.tight_layout()
    outpath = os.path.join(args.outdir, "combined_plots.png")
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved combined plots to {outpath}")


if __name__ == "__main__":
    main()
