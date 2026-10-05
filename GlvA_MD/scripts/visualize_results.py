#!/usr/bin/env python3
"""
Generate a hardcoded PyMOL .pml script from MD pipeline outputs.

This is a plain Python 3 script — no PyMOL installation needed to run it.
It reads the pipeline CSVs, combines them with user-defined residue lists,
and writes a self-contained PyMOL command script.

Usage:
    python3 scripts/visualize_results.py \
        --pdb step5_1.pdb --traj tmp/4_traj_fit.xtc \
        --ligand_resids "450 447 448 449" \
        --internal_water_ids "518 519 528 ..." \
        --user_resids "93 109 110 ..." \
        --csv_ligand_pho results/hbond_PHO.csv \
        --csv_ligand_sugar results/hbond_4GA_0GA_ROH.csv \
        --csv_water results/hbond_water.csv \
        --csv_internal_water results/hbond_internal_water.csv \
        --csv_matrix results/hbond_matrix.csv \
        --output results/visualize_scene.pml

Then open in PyMOL:
    pymol results/visualize_scene.pml
"""

import argparse
import csv
import os
import re
import sys


# ============================================================================
# CSV parsing utilities
# ============================================================================

def _parse_resid_atom(label):
    if "@" not in label:
        return None, None, None
    prefix, atom = label.rsplit("@", 1)
    m = re.match(r"^(\d+)(\D.*)$", prefix)
    if not m:
        return None, None, None
    return m.group(1), m.group(2), atom


def parse_ligand_hbond_csv(csv_path):
    protein_resids = set()
    water_resids = set()

    if not os.path.isfile(csv_path):
        print(f"  [WARN] CSV not found: {csv_path}  — skipping.")
        return protein_resids, water_resids

    with open(csv_path, "r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            for col in ("Donor", "Acceptor"):
                label = row.get(col, "")
                resid, resname, _ = _parse_resid_atom(label)
                if resid is None:
                    continue
                if resname == "WAT":
                    water_resids.add(resid)
                else:
                    protein_resids.add(resid)

    print(f"  Parsed {csv_path}: {len(protein_resids)} protein residues, "
          f"{len(water_resids)} water residues")
    return protein_resids, water_resids


def parse_water_hbond_csv(csv_path):
    water_resids = set()
    protein_resids = set()

    if not os.path.isfile(csv_path):
        print(f"  [WARN] CSV not found: {csv_path}  — skipping.")
        return water_resids, protein_resids

    with open(csv_path, "r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            water_label = row.get("Water ID", "")
            parts = water_label.rsplit("_", 1)
            if len(parts) == 2 and parts[0] == "WAT" and parts[1].isdigit():
                water_resids.add(parts[1])

            prot_label = row.get("Protein Atom", "")
            resid, resname, _ = _parse_resid_atom(prot_label)
            if resid is not None:
                protein_resids.add(resid)

    print(f"  Parsed {csv_path}: {len(water_resids)} unique waters, "
          f"{len(protein_resids)} protein residues")
    return water_resids, protein_resids


def parse_hbond_matrix_csv(csv_path):
    protein_resids = set()

    if not os.path.isfile(csv_path):
        print(f"  [WARN] CSV not found: {csv_path}  — skipping.")
        return protein_resids

    with open(csv_path, "r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            for col in ("hydrogen_label", "acceptor_label"):
                label = row.get(col, "")
                resid, resname, _ = _parse_resid_atom(label)
                if resid is not None:
                    protein_resids.add(resid)

    print(f"  Parsed {csv_path}: {len(protein_resids)} residues")
    return protein_resids


# ============================================================================
# PyMOL selection expression builder
# ============================================================================

def _resids_to_pymol_sel(resids):
    ids = sorted(set(int(r) for r in resids if str(r).lstrip("-").isdigit()))
    if not ids:
        return ""
    return "resi " + "+".join(str(i) for i in ids)


# ============================================================================
# .pml script generation
# ============================================================================

def _emit(w, line=""):
    w.append(line)


def _emit_section(w, title):
    _emit(w)
    _emit(w, f"# ---- {title} ----")


def _emit_cmd(w, cmd_str):
    _emit(w, cmd_str)


def generate_pml(pdb_path, traj_path, ligand_resids, internal_water_resids,
                 csv_water_resids, all_relevant_resids, water_resname,
                 pho_resids=None, lig_resids=None, roh_resids=None,
                 output_png=""):
    
    if pho_resids is None: pho_resids = set()
    if lig_resids is None: lig_resids = set()
    if roh_resids is None: roh_resids = set()
    w = []
    E = lambda line="": _emit(w, line)
    S = lambda title: _emit_section(w, title)
    C = lambda cmd_str: _emit_cmd(w, cmd_str)

    # ---------- Header ----------
    E("# Auto-generated PyMOL scene — GlvA MD analysis pipeline")
    E("from pymol import cmd")
    E("cmd.reinitialize()")
    E()

    # ---------- Load structure and trajectory ----------
    S("Load structure and trajectory")
    C(f'cmd.load("{pdb_path}", "system")')
    if traj_path:
        C(f'cmd.load_traj("{traj_path}", "system")')
    E()

    # ---------- Hide solvent and buffer region ----------
    S("Hide solvent and buffer region")
    C('cmd.hide("everything", "solvent")')
    C('cmd.hide("everything", "resi 451-517")')
    E()

    # ========================================================================
    # Selection groups
    # ========================================================================

    selections_defined = []

    # -- 1. Ligand --
    lig_expr = _resids_to_pymol_sel(ligand_resids)
    if lig_expr:
        S("Selection: ligand")
        C(f'cmd.select("ligand", "{lig_expr}")')
        C('cmd.show("sticks", "ligand")')
        C('cmd.color("green", "ligand")')
        C('cmd.util.cnc("ligand")') # Раскраска гетероатомов
        selections_defined.append("ligand")

        pho_expr = _resids_to_pymol_sel(pho_resids)
        if pho_expr:
            C(f'cmd.select("PHO", "{pho_expr}")')
            C('cmd.show("sticks", "PHO")')
            C('cmd.color("green", "PHO")')
            C('cmd.util.cnc("PHO")')
            selections_defined.append("PHO")

        lig_sub_expr = _resids_to_pymol_sel(lig_resids)
        if lig_sub_expr:
            C(f'cmd.select("LIG", "{lig_sub_expr}")')
            C('cmd.show("sticks", "LIG")')
            C('cmd.color("orange", "LIG")')
            C('cmd.util.cnc("LIG")')
            selections_defined.append("LIG")

        roh_expr = _resids_to_pymol_sel(roh_resids)
        if roh_expr:
            C(f'cmd.select("ROH", "{roh_expr}")')
            C('cmd.show("sticks", "ROH")')
            C('cmd.color("yellow", "ROH")')
            C('cmd.util.cnc("ROH")')
            selections_defined.append("ROH")
        E()

    # -- 2. Internal water --
    iw_expr = _resids_to_pymol_sel(internal_water_resids)
    if iw_expr:
        iw_full = f'resname {water_resname} and ({iw_expr})'
        S("Selection: internal_water")
        C(f'cmd.select("internal_water", "{iw_full}")')
        C('cmd.show("sticks", "internal_water")')
        C('cmd.color("cyan", "internal_water")')
        C('cmd.util.cnc("internal_water")')
        selections_defined.append("internal_water")
        E()

    # -- 3. Unique water --
    pure_unique = csv_water_resids - internal_water_resids
    uw_expr = _resids_to_pymol_sel(pure_unique)
    if uw_expr:
        uw_full = f'resname {water_resname} and ({uw_expr})'
        S("Selection: unique_water")
        C(f'cmd.select("unique_water", "{uw_full}")')
        C('cmd.show("sticks", "unique_water")')
        C('cmd.color("yellow", "unique_water")')
        C('cmd.util.cnc("unique_water")')
        selections_defined.append("unique_water")
        E()

    # -- 4. Relevant residues --
    rel_expr = _resids_to_pymol_sel(all_relevant_resids)
    if rel_expr:
        rel_full = f'not resname {water_resname} and ({rel_expr})'
        S("Selection: relevant_residues")
        C(f'cmd.select("relevant_residues", "{rel_full}")')
        C('cmd.show("sticks", "relevant_residues")')
        C('cmd.color("magenta", "relevant_residues")')
        C('cmd.util.cnc("relevant_residues")')
        # Подписываем ТОЛЬКО релевантные остатки
        C('cmd.label("relevant_residues and name CA", "\'%s %s\' % (resn, resi)")')
        selections_defined.append("relevant_residues")
        E()

    # -- Context --
    S("Context representation")
    # Добавлено 'and not solvent', чтобы скрыть всю остальную воду
    sel_names = " or ".join(selections_defined) if selections_defined else "none"
    C(f'cmd.show("lines", "system and not solvent and not ({sel_names})")')
    C('cmd.color("grey80", "system and not solvent and not ({})")'.format(sel_names))
    E()

    # ========================================================================
    # H-bond analysis
    # ========================================================================

    S("Hydrogen bond analysis")

    C('cmd.set("h_bond_max_angle", 40)')
    C('cmd.set("h_bond_cutoff_center", 3.6)')
    E()

    hbond_targets = [
        ("ligand", "Ligand ↔ Environment"),
        ("internal_water", "Internal Water ↔ Environment"),
        ("unique_water", "Unique Water ↔ Environment"),
    ]

    for target_sel, label in hbond_targets:
        if target_sel in selections_defined:
            env_name = f"env_{target_sel}"
            dist_name = f"hbonds_{target_sel}"

            E(f"# {label}")
            # Если целевая группа - вода, жестко исключаем весь растворитель из ее окружения
            if "water" in target_sel:
                C(f'cmd.select("{env_name}", "byres ({target_sel} around 5) and not solvent")')
            else:
                C(f'cmd.select("{env_name}", "byres ({target_sel} around 5) and not {target_sel}")')
            
            C(f'cmd.distance("{dist_name}", "{target_sel}", "{env_name}", mode=2)')
            E()
    # ========================================================================
    # Final view
    # ========================================================================

    S("Final view")
    
    # Скрываем абсолютно все водороды для чистоты картинки
    C('cmd.hide("sticks", "elem H")')
    C('cmd.hide("lines", "elem H")')

    # Включаем полярные водороды (связанные с N, O или S) только для выделенных групп
    if selections_defined:
        sel_names = " or ".join(selections_defined)
        C(f'cmd.show("sticks", "({sel_names}) and elem H and (neighbor elem N+O+S)")')

    focus = "ligand" if "ligand" in selections_defined else \
            (selections_defined[0] if selections_defined else "system")
    C(f'cmd.zoom("{focus}", 8)')
    C(f'cmd.center("{focus}")')
    C(f'cmd.bg_color("white")')
    C(f'cmd.show("sticks", "resn ROH or PHO or LIG")')
    E()

    if output_png:
        E("# ---- Render PNG (uncomment to enable) ----")
        E(f"# cmd.ray(1920, 1080)")
        E(f'# cmd.png("{output_png}", width=1920, height=1080, dpi=150, ray=1)')
        E()

    C('cmd.orient()')
    return "\n".join(w)


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate a hardcoded PyMOL .pml script from MD pipeline outputs."
    )
    parser.add_argument("--pdb", required=True, help="PDB topology file")
    parser.add_argument("--exclude_resids", default="", help="Остатки, которые нужно исключить из relevant_residues")
    parser.add_argument("--traj", default="", help="Trajectory file (.xtc)")
    parser.add_argument("--ligand_resids", default="", help="Space-separated ligand residue IDs")
    parser.add_argument("--pho_resids", default="", help="Explicit PHO resids")
    parser.add_argument("--lig_resids", default="", help="Explicit LIG resids")
    parser.add_argument("--roh_resids", default="", help="Explicit ROH resids")
    parser.add_argument("--internal_water_ids", default="", help="Internal water resids")
    parser.add_argument("--user_resids", default="", help="User-defined starting resids")
    parser.add_argument("--water_resname", default="WAT", help="Water residue name")
    parser.add_argument("--csv_ligand_pho", default="", help="Path to hbond_PHO.csv")
    parser.add_argument("--csv_ligand_sugar", default="", help="Path to hbond_4GA_0GA_ROH.csv")
    parser.add_argument("--csv_water", default="", help="Path to hbond_water.csv")
    parser.add_argument("--csv_internal_water", default="", help="Path to hbond_internal_water.csv")
    parser.add_argument("--csv_matrix", default="", help="Path to hbond_matrix.csv")
    parser.add_argument("--output", default="results/visualize_scene.pml", help="Output .pml script path")
    parser.add_argument("--png", default="", help="Optional PNG path")
    args = parser.parse_args()

    ligand_resids = set(args.ligand_resids.split()) if args.ligand_resids.strip() else set()
    pho_resids = set(args.pho_resids.split()) if args.pho_resids.strip() else set()
    lig_resids = set(args.lig_resids.split()) if args.lig_resids.strip() else set()
    roh_resids = set(args.roh_resids.split()) if args.roh_resids.strip() else set()
    internal_water_resids = set(args.internal_water_ids.split()) if args.internal_water_ids.strip() else set()
    user_resids = set(args.user_resids.split()) if args.user_resids.strip() else set()
    exclude_resids = set(args.exclude_resids.split()) if args.exclude_resids.strip() else set()

    csv_protein_resids = set()
    csv_water_resids = set()
    csv_internal_protein_resids = set()

    for path in (args.csv_ligand_pho, args.csv_ligand_sugar):
        if path:
            prots, wats = parse_ligand_hbond_csv(path)
            csv_protein_resids |= prots
            csv_water_resids |= wats

    if args.csv_water:
        unique_w, water_ps = parse_water_hbond_csv(args.csv_water)
        csv_water_resids |= unique_w
        csv_protein_resids |= water_ps

    if args.csv_internal_water:
        int_ps, int_ws = parse_ligand_hbond_csv(args.csv_internal_water)
        csv_internal_protein_resids |= int_ps
        csv_water_resids |= int_ws

    if args.csv_matrix:
        mat_ps = parse_hbond_matrix_csv(args.csv_matrix)
        csv_protein_resids |= mat_ps

    all_relevant_resids = (user_resids - exclude_resids | csv_protein_resids | csv_internal_protein_resids)

    pml_content = generate_pml(
        pdb_path=args.pdb,
        traj_path=args.traj,
        ligand_resids=ligand_resids,
        internal_water_resids=internal_water_resids,
        csv_water_resids=csv_water_resids,
        all_relevant_resids=all_relevant_resids,
        water_resname=args.water_resname,
        pho_resids=pho_resids,
        lig_resids=lig_resids,
        roh_resids=roh_resids,
        output_png=args.png,
    )

    output_path = args.output
    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        f.write(pml_content)

if __name__ == "__main__":
    main()
