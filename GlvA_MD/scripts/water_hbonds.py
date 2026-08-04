#!/usr/bin/env python3
"""
Water-mediated H-bond analysis using MDAnalysis.
1) General water analysis: protein residues against all bulk water (WAT).
2) Internal water analysis: protein pocket residues against specific water residue IDs.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import analyze_water_hbonds, analyze_protein_ligand_hbonds

try:
    import MDAnalysis as mda
except ImportError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Water-mediated H-bond analysis"
    )
    parser.add_argument("--tpr", required=True, help="GROMACS topology (.tpr)")
    parser.add_argument("--traj", required=True, help="Trajectory file (.xtc)")
    parser.add_argument("--residues", nargs="+", required=True,
                        help="Protein residue IDs for bulk water analysis")
    parser.add_argument("--water-resname", default="WAT",
                        help="Water residue name (default: WAT)")
    parser.add_argument("--dist-cutoff", type=float, default=3.5,
                        help="Donor-Acceptor distance cutoff")
    parser.add_argument("--angle-cutoff", type=float, default=150,
                        help="D-H-A angle cutoff (degrees)")
    parser.add_argument("--occ-cutoff", type=float, default=5.0,
                        help="Occupancy cutoff (%) for reporting")
    parser.add_argument("--internal-residues", nargs="*", default=[],
                        help="Protein residue IDs for internal water analysis")
    parser.add_argument("--internal-water-ids", nargs="*", default=[],
                        help="Specific water resid IDs for internal water analysis")
    parser.add_argument("--outdir", default="results", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    u = mda.Universe(args.tpr, args.traj)
    print(f"Loaded Universe with {len(u.trajectory)} frames")

    water_resids = [int(r) for r in args.residues]

    print(f"\n{'#'*60}")
    print(f"# 1. Bulk water analysis")
    print(f"{'#'*60}")
    analyze_water_hbonds(
        u, water_resids,
        water_resname=args.water_resname,
        d_a_cutoff=args.dist_cutoff,
        d_h_a_angle_cutoff=args.angle_cutoff,
        cutoff_occ=args.occ_cutoff,
        outdir=args.outdir
    )

    if args.internal_residues and args.internal_water_ids:
        print(f"\n{'#'*60}")
        print(f"# 2. Internal water analysis (specific water molecules)")
        print(f"{'#'*60}")
        internal_resids = [int(r) for r in args.internal_residues]
        internal_waters = [int(w) for w in args.internal_water_ids]
        df_internal = analyze_protein_ligand_hbonds(
            u, internal_resids, internal_waters,
            d_a_cutoff=args.dist_cutoff,
            d_h_a_angle_cutoff=args.angle_cutoff
        )
        csv_path = os.path.join(args.outdir, "hbond_internal_water.csv")
        df_internal.to_csv(csv_path, index=False)
        print(f"Saved to {csv_path}")
    else:
        print("Skipping internal water analysis (no --internal-residues/--internal-water-ids).")


if __name__ == "__main__":
    main()
