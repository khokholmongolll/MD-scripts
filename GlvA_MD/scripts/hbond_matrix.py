#!/usr/bin/env python3
"""
Analyse GROMACS hbond output: parse .xpm + .ndx + .gro,
generate a heatmap and a CSV table of H-bonds above the occupancy cutoff.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import analyze_hbond_matrix


def main():
    parser = argparse.ArgumentParser(description="GROMACS H-bond matrix analysis")
    parser.add_argument("--gro", required=True, help=".gro structure file")
    parser.add_argument("--ndx", required=True, help="hbond.ndx from gmx hbond")
    parser.add_argument("--xpm", required=True, help="hbmat.xpm from gmx hbond")
    parser.add_argument("--cutoff", type=float, default=20.0,
                        help="Occupancy cutoff (%) for filtering")
    parser.add_argument("--outdir", default="results", help="Output directory")
    args = parser.parse_args()

    df = analyze_hbond_matrix(
        gro_file=args.gro,
        ndx_file=args.ndx,
        xpm_file=args.xpm,
        cutoff=args.cutoff,
        outdir=args.outdir
    )

    if df.empty:
        print("No bonds passed the occupancy filter.", file=sys.stderr)


if __name__ == "__main__":
    main()
