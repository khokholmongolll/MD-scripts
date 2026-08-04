# MD Analysis Pipeline — GlvA

Automated, modular Molecular Dynamics analysis pipeline for GROMACS trajectories, driven by a single master bash script with smart incremental caching.

## Directory Structure

```
cryst_wat_opc_ff19SB_5_step/
├── run_analysis.sh              ← Master orchestrator (the only script you run)
├── scripts/
│   ├── utils.py                 ← Shared library: all parsers, analysis functions, plotters
│   ├── hbond_matrix.py          ← GROMACS hbond .xpm → heatmap + occupancy table
│   ├── rmsd_analysis.py         ← RMSD for every index group, overlay + grid plots
│   ├── bond_distance.py         ← Bond distance time-series for specified atom pairs
│   ├── protein_ligand_hbonds.py ← MDAnalysis: protein↔ligand H-bonds + distance profiles
│   ├── water_hbonds.py          ← Bulk water + internal water H-bond analysis
│   ├── combined_plots.py        ← Merged RMSD + bond distance stacked figure
│   └── generate_report.py       ← Assembles all outputs → analysis_report.md
├── tmp/                         ← Intermediate files (.xtc, .ndx, .xpm, .xvg, .state)
└── results/                     ← Final outputs (.png, .csv, analysis_report.md)
```

## Quick Start

```bash
./run_analysis.sh           # normal run (skips cached steps)
./run_analysis.sh --force   # force re-run of every step
```

## Pipeline Overview

The pipeline runs **11 steps** sequentially. Each step's outputs feed into downstream steps.

```
 ┌─────────────────────┐
 │ Step 1  (commented) │  gmx trjconv: PBC centering
 └────────┬────────────┘
          │ traj_pbc.xtc
 ┌────────▼────────────┐
 │ Step 2  (commented) │  gmx trjconv: fit rot+trans → traj_fit.xtc
 └────────┬────────────┘
          │ traj_fit.xtc (= $TRAJ_ANALYSIS)
 ┌────────▼────────────┐
 │ Step 3              │  gmx make_ndx: custom groups → index_updated.ndx
 └────────┬────────────┘
          │ index_updated.ndx
 ┌────────▼────────────┐
 │ Step 4              │  gmx hbond: M6P vs pocket → hbmat.xpm, hbond.ndx
 └────────┬────────────┘
          │ hbmat.xpm, hbond.ndx
 ┌────────▼────────────┐
 │ Step 5              │  Python: parse .xpm, filter by occupancy, draw heatmap
 └───────────┬─────────┘
             │ hbond_matrix.png, hbond_matrix.csv
 ┌────────────▼────────────┐
 │ Step 6                  │  Python (MDAnalysis): RMSD for every index group
 └───────────┬─────────────┘
             │ rmsd_overlay.png, rmsd_subplots.png, rmsd_data.csv
 ┌────────────▼────────────┐
 │ Step 7                  │  Python (MDAnalysis): bond distance time-series
 └───────────┬─────────────┘
             │ bond_distances.png, bond_distances.csv
 ┌────────────▼────────────┐
 │ Step 8                  │  Python (MDAnalysis): Protein↔PHO + Protein↔sugars H-bonds
 └───────────┬─────────────┘
             │ hbond_PHO.csv, hbond_4GA_0GA_ROH.csv, hbond_dist_*.png
 ┌────────────▼────────────┐
 │ Step 9                  │  Python (MDAnalysis): bulk + internal water H-bonds
 └───────────┬─────────────┘
             │ hbond_water.csv, hbond_internal_water.csv
 ┌────────────▼────────────┐
 │ Step 10                 │  Python: stacked figure from RMSD + bond distance CSVs
 └───────────┬─────────────┘
             │ combined_plots.png
 ┌────────────▼────────────┐
 │ Step 11                 │  Python: aggregate all .csv + .png → analysis_report.md
 └─────────────────────────┘
```

## Configuring the Pipeline

All configuration lives at the top of **`run_analysis.sh`**. Edit any variable and re-run — only the affected steps will recalculate.

### Key Variables

| Variable | Default | Description |
|---|---|---|
| `TPR`, `GRO`, `XTC_RAW` | `step5_1.*` | Input topology, coordinates, raw trajectory |
| `NDX_BASE` | `index_new.ndx` | Base GROMACS index file |
| `TRAJ_ANALYSIS` | `$TMP_DIR/4_traj_fit.xtc` | Pre-processed trajectory for analysis |
| `HBOND_GROUP_A`, `HBOND_GROUP_B` | `25`, `36` | Groups for `gmx hbond` |
| `OCCUPANCY_CUTOFF` | `20.0` | Occupancy cutoff (%) for hbond matrix filtering |
| `HBOND_DIST_CUTOFF` | `3.5` | Donor-Acceptor distance (Å) |
| `HBOND_ANGLE_CUTOFF` | `150` | D-H-A angle (degrees) |
| `MIN_OCCUPANCY_BONDS` | `5.0` | Min occupancy (%) for distance profile plots |
| `PROTEIN_RESIDS` | `93 109 ...` | Residue IDs for protein-ligand H-bond analysis |
| `LIGANDS_PHO` | `450` | PHO ligand resname(s) |
| `LIGANDS_SUGARS` | `447 448 449` | Sugar ligand resnames |
| `WATER_PROTEIN_RESIDS` | `93 109 110 170 263 108` | Residues for bulk water analysis |
| `WATER_CUTOFF_OCC` | `10.0` | Occupancy cutoff (%) for bulk water table |
| `INTERNAL_RESIDS` | (large list) | Pocket residues for internal water analysis |
| `INTERNAL_WATER_IDS` | `518 519 ...` | Specific water residue IDs |
| `INTERNAL_WATER_OCC_CUTOFF` | `5.0` | Occupancy cutoff (%) for internal water table |
| `BOND_PAIRS` | `(2647 7039) (4108 7043)` | 0-based atom index pairs |
| `BOND_LABELS` | `ASH_170_OD2...`, `TYD_263_OH...` | Labels for bond distance plots |
| `SKIP_RMSD_GROUPS` | `System Water SOL ...` | Index groups to exclude from RMSD |

### Skipping PBC/Fit Steps

Steps 1 (PBC) and 2 (fit) are commented out by default — the pipeline expects a pre-processed trajectory at `$TRAJ_ANALYSIS`. To enable them, remove the `: <<'SKIP_STEP1'` … `SKIP_STEP1` block delimiters. If you enable them, update `TRAJ_ANALYSIS="$TRAJ_FIT"` to use the freshly generated trajectory.

## Caching System

Every step writes a state fingerprint to `tmp/.step_N.state`. The fingerprint is an md5 hash of:

- **File dependencies** — the md5 checksum of each input file
- **Variable dependencies** — `NAME=value` strings for each config variable

Example state for Step 5:
```
step_hash  step5_1.gro  tmp/hbond.ndx  tmp/hbmat.xpm  "OCC_CUTOFF=20.0"
```

### How Skipping Works

On each run, `skip_step` recomputes the fingerprint and compares with the saved `.state` file. A step is skipped only if **both** of these hold:

1. All expected output files exist
2. The current fingerprint matches the saved one

### Dependency Chain Propagation

Each step's dependencies include the *output files of upstream steps*. If an upstream step re-runs and modifies its output, the downstream step's fingerprint will include a changed file checksum → force it to re-run as well.

**Example**: Change `OCCUPANCY_CUTOFF` from 20.0 to 30.0.
- Step 5's variable dependency `OCC_CUTOFF=30.0` differs from saved `OCC_CUTOFF=20.0` → Step 5 re-runs
- Steps 6–9 are unaffected (they don't depend on `OCCUPANCY_CUTOFF` or Step 5's outputs)
- Step 11 re-runs because its `RESULTS_SNAPSHOT` (a hash of all results/ files) changed

**Example**: Change `HBOND_DIST_CUTOFF` from 3.5 to 3.0.
- Step 8's dependency `DIST_CUT=3.0` differs → Step 8 re-runs
- Step 9's dependency `W_DIST=3.0` differs → Step 9 re-runs
- Steps 3–7 are unaffected
- Step 11 re-runs

Use `--force` to bypass all caching:

```bash
./run_analysis.sh --force
```

## Script Internals

### `scripts/utils.py`

Shared module imported by all Python analysis scripts. Contains:

| Function | Purpose |
|---|---|
| `load_gro_atom_mapping(gro)` | Parse `.gro` → `{atom_num: "ResAtom@Name"}` |
| `parse_gro_mapping(gro)` | Parse `.gro` → `{"ResAtom@Name": 0_based_idx}` |
| `parse_ndx(filename)` | Parse `.ndx` → `{group_name: [0_based_indices]}` |
| `load_hbond_pairs_from_ndx(filename)` | Parse `hbond.ndx` → `{bond_idx: (H_atom, A_atom)}` |
| `parse_xpm(filename)` | Parse GROMACS `.xpm` matrix → binary numpy array |
| `analyze_hbond_matrix(...)` | Full hbond matrix pipeline: parse, filter by occupancy, draw heatmap, save CSV |
| `analyze_protein_ligand_hbonds(u, ...)` | Run two-round MDAnalysis `HydrogenBondAnalysis`, return occupancy table |
| `analyze_water_hbonds(u, ...)` | Bulk water H-bond identification with residue-level tracking |
| `calculate_distances_by_indices(u, ...)` | Compute atom-pair distances across the full trajectory |
| `plot_individual_bonds(times, dists, ...)` | Generate one-subplot-per-bond distance vs. time figure |
| `df_to_bond_list(df, min_occ)` | Convert H-bond DataFrame → list of (Donor, Acceptor) tuples |
| `_to_ns(arr)` | Convert ps → ns for RMSD time axis |

### `scripts/hbond_matrix.py`

**Dependencies**: `--gro`, `--ndx` (hbond.ndx), `--xpm` (hbmat.xpm), `--cutoff`

1. Loads `.gro` atom names, `hbond.ndx` pair mapping, and `.xpm` existence matrix.
2. Filters bonds by occupancy > `--cutoff`.
3. Draws a seaborn heatmap sorted by occupancy.
4. Prints a text report and saves `hbond_matrix.png` + `hbond_matrix.csv`.

### `scripts/rmsd_analysis.py`

**Dependencies**: `--tpr`, `--traj`, `--ndx`, `--skip` (optional)

1. Loads MDAnalysis `Universe`.
2. Parses every group from the `.ndx` file (minus any in `--skip`).
3. Computes RMSD for each group with `superposition_selection="protein"` and `ref_frame=0`.
4. Outputs:
   - `rmsd_overlay.png` — all groups on one figure
   - `rmsd_subplots.png` — 2-column grid of individual subplots
   - `rmsd_data.csv` — `name,time_ps,rmsd`

### `scripts/bond_distance.py`

**Dependencies**: `--tpr`, `--traj`, `--bonds IDX1 IDX2 [...]`, `--labels [...]`

1. Loads MDAnalysis `Universe`.
2. Iterates all frames, computes Euclidean distances for each atom pair.
3. Outputs:
   - `bond_distances.png` — one subplot per bond
   - `bond_distances.csv` — `time_ns,Bond1,Bond2,...`

### `scripts/protein_ligand_hbonds.py`

**Dependencies**: `--tpr`, `--traj`, `--gro`, `--residues`, `--ligands-rounds [...]`, `--dist-cutoff`, `--angle-cutoff`, `--min-occupancy`

1. Loads MDAnalysis `Universe`.
2. For each `--ligands-rounds`, runs `analyze_protein_ligand_hbonds` (two rounds: protein→ligand + ligand→protein).
3. Saves occupancy tables as `hbond_<label>.csv`.
4. For each ligand round, filters bonds above `--min-occupancy`, computes distance profiles via `calculate_distances_by_indices`, and generates individual bond distance plots.

### `scripts/water_hbonds.py`

**Dependencies**: `--tpr`, `--traj`, `--residues`, `--water-resname`, `--occ-cutoff`, `--internal-residues`, `--internal-water-ids`, `--internal-occ-cutoff`

1. **Bulk water analysis** — `analyze_water_hbonds`: identifies water molecules hydrogen-bonded to the specified protein residues above `--occ-cutoff`. Saves `hbond_water.csv`.
2. **Internal water analysis** — `analyze_protein_ligand_hbonds` with specific water residue IDs as "ligands". Filters results to `Occupancy >= --internal-occ-cutoff`. Saves `hbond_internal_water.csv`.

### `scripts/combined_plots.py`

**Dependencies**: (reads data from `results/`)

1. Reads `rmsd_data.csv` and `bond_distances.csv`.
2. Produces a single vertically-stacked figure with all RMSD curves + all bond distance curves.
3. Saves `combined_plots.png`.

### `scripts/generate_report.py`

**Dependencies**: (reads all files from `results/`)

1. Scans `results/` for all `.png` and `.csv` files.
2. Embeds images under section headers (1–5).
3. Converts CSV tables to Markdown tables, rounding `Occupancy` columns to 2 decimal places.
4. Sorts the Bulk Water H-Bonds table by occupancy descending.
5. Skips `rmsd_data.csv`, `bond_distances.csv`, and `combined_plots.png` (these are intermediate data, not report tables).
6. Outputs `analysis_report.md` — open it in any Markdown viewer.

## Output Files

### `results/` (final outputs)

| File | Description |
|---|---|
| `analysis_report.md` | Complete Markdown report with images + tables |
| `hbond_matrix.png` | H-bond existence heatmap |
| `hbond_matrix.csv` | H-bond occupancy table |
| `rmsd_overlay.png` | All-group RMSD overlay |
| `rmsd_subplots.png` | Per-group RMSD grid |
| `rmsd_data.csv` | Raw RMSD values |
| `bond_distances.png` | Bond distance time-series |
| `bond_distances.csv` | Raw bond distance values |
| `hbond_PHO.csv` | PHO ligand H-bond table |
| `hbond_4GA_0GA_ROH.csv` | Sugar ligand H-bond table |
| `hbond_dist_PHO.png` | PHO bond distance dynamics |
| `hbond_dist_4GA_0GA_ROH.png` | Sugar bond distance dynamics |
| `hbond_water.csv` | Bulk water H-bond table |
| `hbond_internal_water.csv` | Internal water H-bond table |
| `combined_plots.png` | Stacked RMSD + bond distances |

### `tmp/` (intermediate, safe to delete)

GROMACS outputs (`.xvg`, `.xpm`, `.ndx`), processed trajectories (`.xtc`), and `.step_N.state` cache files.

## Requirements

- **GROMACS** — `gmx` on `PATH`
- **Python 3** with packages:
  - `MDAnalysis`
  - `numpy`, `pandas`
  - `matplotlib`, `seaborn`
  - `tqdm`

All analysis scripts set `matplotlib.use("Agg")` — they run headless, no display needed.

## Customisation

1. **Change a cutoff**: edit the variable at the top of `run_analysis.sh`, re-run. Only the step that depends on it re-runs.
2. **Add a new residue to the H-bond analysis**: add it to `PROTEIN_RESIDS`, re-run. Step 8 re-runs.
3. **Change bond distance pairs**: edit `BOND_PAIRS` and `BOND_LABELS`, re-run. Step 7 re-runs.
4. **Swap the trajectory**: point `TRAJ_ANALYSIS` to a different `.xtc`, re-run. Steps 4–11 re-run.
5. **Enable PBC/fit**: uncomment Steps 1–2, set `TRAJ_ANALYSIS="$TRAJ_FIT"`. Full pipeline from raw `.xtc`.
