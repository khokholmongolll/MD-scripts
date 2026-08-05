#!/usr/bin/env python3
"""
Aggregate all results into a single Markdown report.
Scans results/ for CSV tables and PNG images, produces analysis_report.md.
"""
import argparse
import os
import glob


def csv_to_markdown_table(csv_path, sort_by=None):
    """Read a CSV and return a Markdown-formatted table string."""
    with open(csv_path, "r") as f:
        lines = [ln.rstrip("\n") for ln in f.readlines()]
    if not lines:
        return "*Empty table.*\n"

    rows = [ln.split(",") for ln in lines]
    header = rows[0]
    data = rows[1:]

    occ_idx = None
    for i, h in enumerate(header):
        if "occupancy" in h.lower():
            occ_idx = i
            break

    if occ_idx is not None:
        for row in data:
            try:
                row[occ_idx] = f"{float(row[occ_idx]):.2f}"
            except (ValueError, IndexError):
                pass

    if sort_by is not None and sort_by in header:
        col = header.index(sort_by)

        def _sort_key(r):
            try:
                return float(r[col])
            except (ValueError, IndexError):
                return float("-inf")

        data.sort(key=_sort_key, reverse=True)

    n_cols = len(header)
    col_widths = [len(h) for h in header]
    for row in data:
        for i in range(min(n_cols, len(row))):
            col_widths[i] = max(col_widths[i], len(str(row[i])))

    header_line = "| " + " | ".join(
        str(h).ljust(col_widths[i]) for i, h in enumerate(header)
    ) + " |"
    sep_line = "|-" + "-|-".join(
        "-" * col_widths[i] for i in range(n_cols)
    ) + "-|"

    out = [header_line, sep_line]
    for row in data:
        cells = [str(row[i]) if i < len(row) else "" for i in range(n_cols)]
        out.append("| " + " | ".join(
            c.ljust(col_widths[i]) for i, c in enumerate(cells)
        ) + " |")
    return "\n".join(out) + "\n"


def main():
    parser = argparse.ArgumentParser(description="Generate final Markdown report")
    parser.add_argument("--outdir", default="results", help="Results directory")
    args = parser.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    report_path = os.path.join(outdir, "analysis_report.md")

    pngs = sorted(glob.glob(os.path.join(outdir, "*.png")))
    csvs = sorted(glob.glob(os.path.join(outdir, "*.csv")))

    sections = []

    sections.append("# MD Analysis Report\n")

    image_section_map = {
        "hbond_matrix.png": "## 1. H-Bond Matrix (GROMACS)",
        "rmsd_overlay.png": "## 2. RMSD Analysis — All Groups",
        "rmsd_subplots.png": "## 3. RMSD — Individual Groups",
        "bond_distances.png": "## 4. Bond Distance Analysis",
    }
    image_skip = {"combined_plots.png"}
    hbond_dist_images = [p for p in pngs if "hbond_dist_" in os.path.basename(p)]
    if hbond_dist_images:
        image_section_map.update({
            os.path.basename(p): "## 5. H-Bond Distance Dynamics"
            for p in sorted(hbond_dist_images)
        })

    # Metadynamics images — placed after the MD sections
    meta_section_base = 6
    image_section_map.update({
        "meta_colvar.png": f"## {meta_section_base}. Metadynamics — CV Trajectories",
        "meta_fes2d.png": f"## {meta_section_base + 1}. Metadynamics — 2D Free Energy Surface",
        "meta_convergence.png": f"## {meta_section_base + 2}. Metadynamics — FES Convergence",
    })

    csv_section_names = {
        "hbond_matrix.csv": "### H-Bond Matrix Table",
        "hbond_PHO.csv": "### PHO H-Bonds",
        "hbond_4GA_0GA_ROH.csv": "### 4GA / 0GA / ROH H-Bonds",
        "hbond_water.csv": "### Bulk Water H-Bonds",
        "hbond_internal_water.csv": "### Internal Water H-Bonds",
        "meta_fes2d.csv": "### 2D FES Grid Data",
    }
    csv_skip = {"rmsd_data.csv", "bond_distances.csv", "meta_fes2d.csv"}

    # ---- Metadynamics skipped? ----
    meta_skipped = os.path.isfile(os.path.join(outdir, ".meta_skipped"))

    seen_sections = set()
    for png in pngs:
        name = os.path.basename(png)
        if name in image_skip:
            continue
        section_title = image_section_map.get(name)
        if section_title and section_title not in seen_sections:
            sections.append(f"\n{section_title}\n")
            seen_sections.add(section_title)
        sections.append(f"![{name}]({name})\n")

    for csv in csvs:
        name = os.path.basename(csv)
        if name in csv_skip:
            continue
        section_title = csv_section_names.get(name, f"### {name}")
        if section_title and section_title not in seen_sections:
            sections.append(f"\n{section_title}\n")
            seen_sections.add(section_title)

        try:
            if name == "hbond_water.csv":
                table = csv_to_markdown_table(csv, sort_by="Occupancy (%)")
            else:
                table = csv_to_markdown_table(csv)
            sections.append(table)
        except Exception as e:
            sections.append(f"*Could not parse {csv}: {e}*\n")

    # Report that metadynamics was skipped when HILLS / COLVAR are absent
    if meta_skipped:
        sections.append("\n## Metadynamics\n")
        sections.append(
            "*Stage skipped: HILLS and/or COLVAR file "
            "not found in the working directory.*\n"
        )

    with open(report_path, "w") as f:
        f.write("\n".join(sections))

    print(f"Report written to {report_path}")
    print(f"  Images embedded: {len(pngs)}")
    print(f"  Tables embedded: {len(csvs)}")


if __name__ == "__main__":
    main()
