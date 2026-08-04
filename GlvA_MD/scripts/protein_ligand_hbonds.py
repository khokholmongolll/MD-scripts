#!/usr/bin/env python3
"""
Protein-ligand H-bond analysis using MDAnalysis.
Runs H-bond detection, generates tables, bond distance plots.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import (analyze_protein_ligand_hbonds, parse_gro_mapping,
                   calculate_distances_by_indices, plot_individual_bonds,
                   df_to_bond_list)

try:
    import MDAnalysis as mda
except ImportError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Protein-ligand H-bond analysis with distance plots"
    )
    parser.add_argument("--tpr", required=True, help="GROMACS topology (.tpr)")
    parser.add_argument("--traj", required=True, help="Trajectory file (.xtc)")
    parser.add_argument("--gro", required=True, help="Coordinate file (.gro)")
    parser.add_argument("--residues", nargs="+", required=True,
                        help="Protein residue IDs")
    parser.add_argument("--ligands-rounds", nargs="+", action="append",
                        help="Ligand resnames for each round: "
                             "--ligands-rounds PHO --ligands-rounds 4GA 0GA ROH")
    parser.add_argument("--round-labels", nargs="*", default=[],
                        help="Labels for each ligand round (for output naming)")
    parser.add_argument("--dist-cutoff", type=float, default=3.5,
                        help="Donor-Acceptor distance cutoff")
    parser.add_argument("--angle-cutoff", type=float, default=150,
                        help="D-H-A angle cutoff (degrees)")
    parser.add_argument("--min-occupancy", type=float, default=5.0,
                        help="Min occupancy (%) for bond distance plots")
    parser.add_argument("--outdir", default="results", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if not args.ligands_rounds:
        print("ERROR: At least one --ligands-rounds required", file=sys.stderr)
        sys.exit(1)

    u = mda.Universe(args.tpr, args.traj)
    print(f"Loaded Universe with {len(u.trajectory)} frames")

    protein_resids = [int(r) for r in args.residues]

    round_info = list(args.ligands_rounds)
    labels = list(args.round_labels)
    while len(labels) < len(round_info):
        labels.append("_".join(round_info[len(labels)]))

    all_bonds = {}
    for ligs, label in zip(round_info, labels):
        print(f"\n{'#'*60}")
        print(f"# Round: {label}")
        print(f"# Ligands: {ligs}")
        print(f"{'#'*60}")
        df = analyze_protein_ligand_hbonds(
            u, protein_resids, ligs,
            d_a_cutoff=args.dist_cutoff,
            d_h_a_angle_cutoff=args.angle_cutoff
        )
        csv_path = os.path.join(args.outdir, f"hbond_{label}.csv")
        df.to_csv(csv_path, index=False)
        print(f"Saved table to {csv_path}")
        all_bonds[label] = df

    gro_mapping = parse_gro_mapping(args.gro)

    for label, df in all_bonds.items():
        bonds = df_to_bond_list(df, min_occupancy=args.min_occupancy)
        if not bonds:
            print(f"No bonds above {args.min_occupancy}% occupancy for '{label}'. "
                  f"Skipping distance plot.")
            continue
        print(f"\nComputing distance profiles for {len(bonds)} bonds in '{label}'...")
        times, dists = calculate_distances_by_indices(u, gro_mapping, bonds)
        if times is not None and dists is not None:
            plot_individual_bonds(times, dists,
                                  title_prefix=f"{label}", outdir=args.outdir)


if __name__ == "__main__":
    main()
