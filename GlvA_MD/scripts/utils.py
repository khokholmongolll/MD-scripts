import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os
from tqdm import tqdm


# ====================================================================
# .gro file parsers
# ====================================================================

def load_gro_atom_mapping(gro_filename):
    atom_map = {}
    print(f"Loading structure from {gro_filename}...")
    with open(gro_filename, 'r') as f:
        lines = f.readlines()
    for line in lines[2:-1]:
        if len(line) < 20:
            continue
        try:
            res_num = line[0:5].strip()
            res_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_num = int(line[15:20].strip())
            atom_map[atom_num] = f"{res_num}{res_name}@{atom_name}"
        except ValueError:
            continue
    print(f"  Loaded {len(atom_map)} atoms from {gro_filename}")
    return atom_map


def parse_gro_mapping(gro_path="step5_1.gro"):
    gro_mapping = {}
    atom_idx_0_based = 0
    with open(gro_path, "r") as f:
        lines = f.readlines()
    for line in lines[2:-1]:
        if len(line.strip()) < 20:
            continue
        res_num = line[0:5].strip()
        res_name = line[5:10].strip()
        atom_name = line[10:15].strip()
        key = f"{res_num}{res_name}@{atom_name}"
        gro_mapping[key] = atom_idx_0_based
        atom_idx_0_based += 1
    print(f"  [GRO] Indexed {len(gro_mapping)} atoms from {gro_path}")
    return gro_mapping


# ====================================================================
# .ndx file parsers
# ====================================================================

def parse_ndx(filename):
    groups = {}
    current_group = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('['):
                current_group = line.strip('[] ').strip()
                groups[current_group] = []
            elif current_group and line:
                groups[current_group].extend([int(x) - 1 for x in line.split()])
    return groups


def load_hbond_pairs_from_ndx(ndx_filename):
    pairs = {}
    current_idx = 0
    in_target_section = False
    with open(ndx_filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('['):
                in_target_section = 'hbonds' in line.lower()
                continue
            if in_target_section:
                parts = line.split()
                if len(parts) >= 3:
                    pairs[current_idx] = (int(parts[1]), int(parts[2]))
                    current_idx += 1
    return pairs


# ====================================================================
# XPM parser
# ====================================================================

def parse_xpm(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    data_lines = [line.split('"')[1] for line in lines if line.startswith('"') and len(line.split('"')) > 1]
    meta = data_lines[0].split()
    y_dim = int(meta[1])
    num_colors = int(meta[2])
    matrix_lines = data_lines[1 + num_colors:]
    binary_matrix = np.zeros((y_dim, len(matrix_lines[0])), dtype=int)
    for i, line in enumerate(matrix_lines):
        for j, char in enumerate(line):
            if char != ' ':
                binary_matrix[i, j] = 1
    return binary_matrix


# ====================================================================
# GROMACS H-bond matrix analysis
# ====================================================================

def analyze_hbond_matrix(gro_file="step5_1.gro", ndx_file="hbond.ndx",
                         xpm_file="hbmat.xpm", cutoff=20.0, outdir="results"):
    os.makedirs(outdir, exist_ok=True)

    gro_atoms = load_gro_atom_mapping(gro_file)
    atom_mapping = load_hbond_pairs_from_ndx(ndx_file)
    hb_matrix_full = parse_xpm(xpm_file)

    occupancy = np.mean(hb_matrix_full, axis=1) * 100
    print(f"Loaded data from {xpm_file}. Filtering at occupancy > {cutoff}%...")

    valid_indices = [i for i, occ in enumerate(occupancy) if occ > cutoff]
    if len(valid_indices) == 0:
        print(f"No bonds exceed {cutoff}% occupancy.")
        return pd.DataFrame()

    hb_matrix_filtered = hb_matrix_full[valid_indices, :]

    report_data = []
    yticklabels_labels = []

    for idx, i in enumerate(valid_indices):
        occ = occupancy[i]
        h_atom_idx, a_atom_idx = atom_mapping.get(i, (None, None))
        h_label = gro_atoms.get(h_atom_idx, f"Atom_{h_atom_idx}")
        a_label = gro_atoms.get(a_atom_idx, f"Atom_{a_atom_idx}")
        report_data.append({
            "bond_idx": i,
            "occupancy": occ,
            "hydrogen_label": h_label,
            "acceptor_label": a_label
        })
        yticklabels_labels.append(f"Idx {i} ({h_label} ... {a_label})")

    df_report = pd.DataFrame(report_data)

    sort_order = np.argsort(df_report["occupancy"])[::-1]
    df_report = df_report.iloc[sort_order].reset_index(drop=True)
    hb_matrix_filtered = hb_matrix_filtered[sort_order, :]
    yticklabels_labels = [yticklabels_labels[j] for j in sort_order]

    plt.figure(figsize=(14, 2 + len(valid_indices) * 0.4))
    ax = plt.gca()
    sns.heatmap(hb_matrix_filtered, cmap="Blues", cbar=False,
                yticklabels=yticklabels_labels, ax=ax)
    plt.title(f"H-Bond Existence Map (Occupancy > {cutoff}%)")
    plt.xlabel("Frame")
    plt.ylabel("Hydrogen Bond (H ... A)")
    plt.yticks(rotation=0)
    plt.tight_layout()
    outpath = os.path.join(outdir, "hbond_matrix.png")
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved heatmap to {outpath}")

    sep = "=" * 75
    print(f"\n{sep}")
    print(f"{f'FINAL HYDROGEN BONDS REPORT (Occupancy > {cutoff}%)':^75}")
    print(sep)
    print(f"{'Index':<6} | {'Occupancy':<9} | {'Hydrogen Group':<20} <...> {'Acceptor Group':<20}")
    print("-" * 75)
    for _, row in df_report.iterrows():
        print(f"{int(row['bond_idx']):<6} | {row['occupancy']:>7.1f}% | "
              f"{row['hydrogen_label']:<20} <...> {row['acceptor_label']:<20}")
    print(sep)

    csv_path = os.path.join(outdir, "hbond_matrix.csv")
    df_report.to_csv(csv_path, index=False)
    return df_report


# ====================================================================
# MDAnalysis H-bond analysis (protein <-> ligand)
# ====================================================================

def analyze_protein_ligand_hbonds(u, protein_resids, ligand_resnames,
                                  d_a_cutoff=3.5, d_h_a_angle_cutoff=150):
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

    print(f"\n=== H-bond analysis ===")
    print(f"Protein residues: {protein_resids}")
    print(f"Ligands: {ligand_resnames}\n")

    protein_query = " or ".join([f"resid {r}" for r in protein_resids])
    ligand_query = " or ".join([f"resid {l}" for l in ligand_resnames])

    hydrogens_prot = f"({protein_query}) and name H*"
    acceptors_lig = f"({ligand_query}) and (name O* or name N*)"
    hydrogens_lig = f"({ligand_query}) and name H*"
    acceptors_prot = f"({protein_query}) and (name O* or name N*)"
    donors_prot = f"({protein_query}) and (name O* or name N*)"
    donors_lig = f"({ligand_query}) and (name O* or name N*)"

    print("1. Round A (Protein -> Ligand)...")
    hb_round_A = HydrogenBondAnalysis(
        universe=u,
        donors_sel=donors_prot,
        hydrogens_sel=hydrogens_prot,
        acceptors_sel=acceptors_lig,
        d_a_cutoff=d_a_cutoff,
        d_h_a_angle_cutoff=d_h_a_angle_cutoff,
        update_selections=False
    )
    hb_round_A.run(verbose=True)

    print("2. Round B (Ligand -> Protein)...")
    hb_round_B = HydrogenBondAnalysis(
        universe=u,
        donors_sel=donors_lig,
        hydrogens_sel=hydrogens_lig,
        acceptors_sel=acceptors_prot,
        d_a_cutoff=d_a_cutoff,
        d_h_a_angle_cutoff=d_h_a_angle_cutoff,
        update_selections=False
    )
    hb_round_B.run(verbose=True)

    print("3. Collecting and sorting results...")
    total_frames = len(u.trajectory)
    bond_counts = {}

    for hb_instance in [hb_round_A, hb_round_B]:
        if hb_instance.results.hbonds is not None and len(hb_instance.results.hbonds) > 0:
            for bond in hb_instance.results.hbonds:
                frame, d_idx, h_idx, a_idx, dist, ang = bond
                pair = (int(d_idx), int(a_idx))
                if pair not in bond_counts:
                    bond_counts[pair] = 0
                bond_counts[pair] += 1

    report = []
    for (d_idx, a_idx), freq in bond_counts.items():
        donor = u.atoms[d_idx]
        acceptor = u.atoms[a_idx]
        occ = (freq / total_frames) * 100
        report.append({
            'Donor': f"{donor.resid}{donor.resname}@{donor.name}",
            'Acceptor': f"{acceptor.resid}{acceptor.resname}@{acceptor.name}",
            'Frequency': freq,
            'Occupancy (%)': occ
        })

    df_hb = pd.DataFrame(report)

    if not df_hb.empty:
        df_hb = df_hb.sort_values(by='Occupancy (%)', ascending=False).reset_index(drop=True)

        sep = "=" * 85
        print(f"\n{sep}")
        print(f"{'ARRAY HYDROGEN BONDS REPORT':^85}")
        print(sep)
        print(f"{'Occupancy':<9} | {'Donor':<30} ---> {'Acceptor':<30}")
        print("-" * 85)
        for _, row in df_hb.iterrows():
            print(f"{row['Occupancy (%)']:>7.1f}% | {row['Donor']:<30} ---> {row['Acceptor']:<30}")
        print(sep + "\n")
    else:
        print("No hydrogen bonds found.\n")

    return df_hb


# ====================================================================
# Water-mediated H-bond analysis
# ====================================================================

def analyze_water_hbonds(u, protein_resids, water_resname="WAT",
                         d_a_cutoff=3.5, d_h_a_angle_cutoff=150,
                         cutoff_occ=5.0, outdir="results"):
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

    os.makedirs(outdir, exist_ok=True)

    print(f"\n=== Water-mediated H-bond analysis ===")
    print(f"Protein residues: {protein_resids}")
    print(f"Water resname: {water_resname}")
    print(f"Occupancy cutoff: {cutoff_occ}%\n")

    protein_query = " or ".join([f"resid {r}" for r in protein_resids])

    donors_prot = f"({protein_query}) and (name O* or name N*)"
    hydrogens_prot = f"({protein_query}) and name H*"
    acceptors_prot = f"({protein_query}) and (name O* or name N*)"

    donors_wat = f"resname {water_resname} and (name O* or name N*)"
    hydrogens_wat = f"resname {water_resname} and name H*"
    acceptors_wat = f"resname {water_resname} and (name O* or name N*)"

    print("1. Round A (Protein -> Water)...")
    hb_round_A = HydrogenBondAnalysis(
        universe=u, donors_sel=donors_prot, hydrogens_sel=hydrogens_prot,
        acceptors_sel=acceptors_wat,
        d_a_cutoff=d_a_cutoff, d_h_a_angle_cutoff=d_h_a_angle_cutoff,
        update_selections=False
    )
    hb_round_A.run(verbose=True)

    print("\n2. Round B (Water -> Protein)...")
    hb_round_B = HydrogenBondAnalysis(
        universe=u, donors_sel=donors_wat, hydrogens_sel=hydrogens_wat,
        acceptors_sel=acceptors_prot,
        d_a_cutoff=d_a_cutoff, d_h_a_angle_cutoff=d_h_a_angle_cutoff,
        update_selections=False
    )
    hb_round_B.run(verbose=True)

    print("\n3. Identifying specific water molecules...")
    total_frames = len(u.trajectory)
    bond_counts = {}

    if hb_round_A.results.hbonds is not None and len(hb_round_A.results.hbonds) > 0:
        for bond in tqdm(hb_round_A.results.hbonds, desc="Round A"):
            frame, d_idx, h_idx, a_idx, dist, ang = bond
            p_idx = int(d_idx)
            w_idx = int(a_idx)
            w_resid = u.atoms[w_idx].resid
            key = (p_idx, w_resid, 'Donor (gives H to water)')
            if key not in bond_counts:
                bond_counts[key] = set()
            bond_counts[key].add(int(frame))

    if hb_round_B.results.hbonds is not None and len(hb_round_B.results.hbonds) > 0:
        for bond in tqdm(hb_round_B.results.hbonds, desc="Round B"):
            frame, d_idx, h_idx, a_idx, dist, ang = bond
            w_idx = int(d_idx)
            p_idx = int(a_idx)
            w_resid = u.atoms[w_idx].resid
            key = (p_idx, w_resid, 'Acceptor (receives H from water)')
            if key not in bond_counts:
                bond_counts[key] = set()
            bond_counts[key].add(int(frame))

    print("Generating final table...")
    report = []
    for (p_idx, w_resid, role), frames in bond_counts.items():
        occ = (len(frames) / total_frames) * 100
        if occ >= cutoff_occ:
            atom = u.atoms[p_idx]
            report.append({
                'Occupancy (%)': occ,
                'Protein Atom': f"{atom.resid}{atom.resname}@{atom.name}",
                'Water ID': f"{water_resname}_{w_resid}",
                'Role': role
            })

    df_hb = pd.DataFrame(report)

    if not df_hb.empty:
        df_hb = df_hb.sort_values(by=['Protein Atom', 'Occupancy (%)'],
                                  ascending=[True, False]).reset_index(drop=True)

        sep = "=" * 85
        print(f"\n{sep}")
        print(f"{f'ORDERED WATER MOLECULES (Occupancy > {cutoff_occ}%)':^85}")
        print(sep)
        print(f"{'Occupancy':<9} | {'Protein Atom':<20} | {'Water ID':<15} | {'Role':<30}")
        print("-" * 85)
        for _, row in df_hb.iterrows():
            print(f"{row['Occupancy (%)']:>7.1f}% | {row['Protein Atom']:<20} | "
                  f"{row['Water ID']:<15} | {row['Role']:<30}")
        print(sep)

        csv_path = os.path.join(outdir, "hbond_water.csv")
        df_hb.to_csv(csv_path, index=False)
        print(f"  Saved to {csv_path}")
    else:
        print(f"No ordered water molecules with occupancy > {cutoff_occ}% found.")

    return df_hb


# ====================================================================
# Trajectory distance calculations
# ====================================================================

def calculate_distances_by_indices(u, gro_mapping, bonds_list):
    valid_pairs = []
    for d_label, a_label in bonds_list:
        if d_label not in gro_mapping or a_label not in gro_mapping:
            print(f"WARNING: Pair {d_label} -> {a_label} not found in .gro. Skipping.")
            continue
        valid_pairs.append(
            (gro_mapping[d_label], gro_mapping[a_label], f"{d_label} ... {a_label}")
        )

    if not valid_pairs:
        print("No valid pairs for distance calculation!")
        return None, None

    distances_dict = {label: [] for _, _, label in valid_pairs}
    times = []

    print(f"  [MD] Computing distances for {len(valid_pairs)} bonds...")
    for ts in tqdm(u.trajectory, desc="Trajectory"):
        times.append(ts.time)
        for idx1, idx2, label in valid_pairs:
            p1 = u.atoms[idx1].position
            p2 = u.atoms[idx2].position
            distances_dict[label].append(np.linalg.norm(p1 - p2))

    times = np.array(times) / 1000.0  # ps -> ns
    return times, distances_dict


def plot_individual_bonds(times, distances_dict, title_prefix="Cluster",
                          outdir="results"):
    os.makedirs(outdir, exist_ok=True)

    labels = list(distances_dict.keys())
    n_plots = len(labels)

    fig, axes = plt.subplots(n_plots, 1, figsize=(14, 3.5 * n_plots), sharex=True)
    if n_plots == 1:
        axes = [axes]

    fig.suptitle(f"Individual Hydrogen Bond Dynamics: {title_prefix}",
                 fontsize=16, fontweight='bold', y=1.01)

    for i, label in enumerate(labels):
        ax = axes[i]
        raw_dist = np.array(distances_dict[label])
        smoothed_dist = pd.Series(raw_dist).rolling(window=10, min_periods=1).mean()

        ax.plot(times, raw_dist, color='gray', alpha=0.25, linewidth=0.8, label="Raw data")
        ax.plot(times, smoothed_dist, color='tab:blue', alpha=0.9, linewidth=1.5, label="Smoothed trend")
        ax.axhline(y=3.5, color='tab:red', linestyle='--', linewidth=1.5, label="Cutoff (3.5 A)")
        ax.set_title(label, fontsize=11, fontweight='semibold', loc='left', color='dimgray')
        ax.set_ylabel("Distance (A)", fontsize=10)
        ax.set_ylim(2.0, 5.5)
        ax.grid(True, linestyle=':', alpha=0.5)
        if i == 0:
            ax.legend(loc="upper right", fontsize=9, frameon=True)

    axes[-1].set_xlabel("Time (ns)", fontsize=12)
    plt.tight_layout()

    safe_name = title_prefix.replace(" & ", "_").replace(" ", "_").replace(",", "")
    outpath = os.path.join(outdir, f"hbond_dist_{safe_name}.png")
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved plot to {outpath}")


def df_to_bond_list(df, min_occupancy=0.0):
    if df is None or df.empty:
        return []
    filtered_df = df[df['Occupancy (%)'] >= min_occupancy]
    return list(zip(filtered_df['Donor'], filtered_df['Acceptor']))


# ====================================================================
# Helpers
# ====================================================================

def _to_ns(arr):
    return arr[:, 1] / 1000.0
