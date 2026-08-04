#!/bin/bash
set -e
# ============================================================
# Master Analysis Pipeline for MD Simulation (GlvA)
# ============================================================

# --------------- Input files ---------------
TPR="step5_1.tpr"
GRO="step5_1.gro"
XTC_RAW="step5_1.xtc"
NDX_BASE="index_new.ndx"

# --------------- Directories ---------------
SCRIPTS_DIR="scripts"
TMP_DIR="tmp"
RESULTS_DIR="results"

# --------------- Trajectory processing ---------------
# Set TRAJ_ANALYSIS to a pre-processed .xtc to skip PBC/fit steps.
TRAJ_PBC="$TMP_DIR/traj_pbc.xtc"
TRAJ_FIT="$TMP_DIR/4_traj_fit.xtc"
TRAJ_ANALYSIS="$TRAJ_FIT"

# --------------- Output files ---------------
NDX_UPDATED="$TMP_DIR/index_updated.ndx"

# --------------- GROMACS groups ---------------
CENTER_GROUP=4        # Backbone
OUTPUT_GROUP=0         # System
FIT_GROUP=4            # Backbone

# --------------- H-bond matrix (gmx hbond) ---------------
HBOND_GROUP_A=25       # M6P in index_updated.ndx
HBOND_GROUP_B=36       # pocket  in index_updated.ndx

# --------------- Analysis parameters ---------------
OCCUPANCY_CUTOFF=20.0
HBOND_DIST_CUTOFF=3.5
HBOND_ANGLE_CUTOFF=150
MIN_OCCUPANCY_BONDS=5.0

# --------------- Protein residues for PDAnalysis H-bond ---------------
PROTEIN_RESIDS="93 109 110 170 263 171 86 262 147 445 261 318 108"
LIGANDS_PHO="450"
LIGANDS_SUGARS="447 448 449"

# --------------- Water analysis ---------------
WATER_RESNAME="WAT"
WATER_CUTOFF_OCC=5.0
WATER_PROTEIN_RESIDS="93 109 110 170 263 108"

# Internal water analysis (specific water residue IDs)
INTERNAL_RESIDS="15 86 108 109 110 111 118 144 145 146 147 148 149 150 151 152 167 168 169 170 171 172 173 174 175 199 200 201 261 262 263 266 318 321 338 340 357"
INTERNAL_WATER_IDS="518 519 528 530 531 535 536 539 541 542 547 550 552 556 558 564 566 575"

# --------------- Bond distance atom pairs (0-based) ---------------
BOND_PAIRS=("2647 7039" "4108 7043")
BOND_LABELS=("ASH_170_OD2--4GA_448_O4" "TYD_263_OH--0GA_449_H2")

# --------------- RMSD skip groups ---------------
# Leave empty to compute all. Example: "System Water SOL non-Protein Water_and_ions K+ Cl-"
SKIP_RMSD_GROUPS="System Water SOL non-Protein Water_and_ions K+ Cl-"

# ============================================================
# Setup
# ============================================================
mkdir -p "$TMP_DIR" "$RESULTS_DIR"

# ============================================================
# Step 1: PBC centering
# ============================================================
#echo "=== Step 1: PBC centering ==="
#echo "$CENTER_GROUP $OUTPUT_GROUP" | gmx trjconv \ 
#    -s "$TPR" -f "$XTC_RAW" -n "$NDX_BASE" \ 
#    -pbc mol -center -ur compact \ 
#    -o "$TRAJ_PBC"

# ============================================================
# Step 2: Fit rotation + translation
# ============================================================
#echo "=== Step 2: Fit rotation + translation ==="
#echo "$FIT_GROUP $OUTPUT_GROUP" | gmx trjconv \
#    -s "$TPR" -f "$TRAJ_PBC" -n "$NDX_BASE" \
#    -fit rot+trans \
#    -o "$TRAJ_FIT"

# ============================================================
# Step 3: Build custom index file
# ============================================================
echo "=== Step 3: Building custom index ==="
gmx make_ndx -f "$GRO" -n "$NDX_BASE" -o "$NDX_UPDATED" > "$TMP_DIR/make_ndx.log" 2>&1 <<'EOF'
1 | 14
r 170
r 263
r 109
r 110
18 | 19
16 | 17
r 14 15 22 86 90 93 97 108 109 110 111 115 118 145 146 147 148 149 168 169 170 171 172 173 174 175 177 178 195 199 200 201 202 204 228 238 242 243 244 245 246 247 248 260 261 262 263 264 265 266 273 276 277 281 283 284 286 287 291 314 315 316 317 318 319 321 357 445
name 36 pocket
r 93
27 | 29
name 38 complex
q
EOF
echo "  Created $NDX_UPDATED"

# ============================================================
# Step 4: GROMACS H-bond analysis (matrix)
# ============================================================
echo "=== Step 4: GROMACS H-bond analysis ==="
echo "$HBOND_GROUP_A $HBOND_GROUP_B" | gmx hbond \
    -f "$TRAJ_ANALYSIS" \
    -s "$TPR" \
    -n "$NDX_UPDATED" \
    -num "$TMP_DIR/hbnum.xvg" \
    -dist "$TMP_DIR/hbdist.xvg" \
    -hbm "$TMP_DIR/hbmat.xpm" \
    -hbn "$TMP_DIR/hbond.ndx" \
    -nomerge

# ============================================================
# Step 5: H-bond matrix analysis + heatmap
# ============================================================
echo "=== Step 5: H-bond matrix analysis ==="
python3 "$SCRIPTS_DIR/hbond_matrix.py" \
    --gro "$GRO" \
    --ndx "$TMP_DIR/hbond.ndx" \
    --xpm "$TMP_DIR/hbmat.xpm" \
    --cutoff "$OCCUPANCY_CUTOFF" \
    --outdir "$RESULTS_DIR"

# ============================================================
# Step 6: RMSD analysis for all index groups
# ============================================================
echo "=== Step 6: RMSD analysis ==="
RMSD_SKIP_ARGS=""
[ -n "$SKIP_RMSD_GROUPS" ] && RMSD_SKIP_ARGS="--skip $SKIP_RMSD_GROUPS"
python3 "$SCRIPTS_DIR/rmsd_analysis.py" \
    --tpr "$TPR" \
    --traj "$TRAJ_ANALYSIS" \
    --ndx "$NDX_UPDATED" \
    $RMSD_SKIP_ARGS \
    --outdir "$RESULTS_DIR"

# ============================================================
# Step 7: Bond distance analysis
# ============================================================
echo "=== Step 7: Bond distance analysis ==="
BOND_ARGS=""
for pair in "${BOND_PAIRS[@]}"; do
    BOND_ARGS="$BOND_ARGS --bonds $pair"
done
LABEL_ARGS=""
for lbl in "${BOND_LABELS[@]}"; do
    LABEL_ARGS="$LABEL_ARGS --labels $lbl"
done
python3 "$SCRIPTS_DIR/bond_distance.py" \
    --tpr "$TPR" \
    --traj "$TRAJ_ANALYSIS" \
    $BOND_ARGS \
    $LABEL_ARGS \
    --outdir "$RESULTS_DIR"

# ============================================================
# Step 8: Protein-ligand H-bond analysis
# ============================================================
echo "=== Step 8: Protein-ligand H-bond analysis ==="
python3 "$SCRIPTS_DIR/protein_ligand_hbonds.py" \
    --tpr "$TPR" \
    --traj "$TRAJ_ANALYSIS" \
    --gro "$GRO" \
    --residues $PROTEIN_RESIDS \
    --ligands-rounds "$LIGANDS_PHO" \
    --ligands-rounds $LIGANDS_SUGARS \
    --round-labels "PHO" "4GA_0GA_ROH" \
    --dist-cutoff "$HBOND_DIST_CUTOFF" \
    --angle-cutoff "$HBOND_ANGLE_CUTOFF" \
    --min-occupancy "$MIN_OCCUPANCY_BONDS" \
    --outdir "$RESULTS_DIR"

# ============================================================
# Step 9: Water-mediated H-bond analysis
# ============================================================
echo "=== Step 9: Water-mediated H-bond analysis ==="
python3 "$SCRIPTS_DIR/water_hbonds.py" \
    --tpr "$TPR" \
    --traj "$TRAJ_ANALYSIS" \
    --residues $WATER_PROTEIN_RESIDS \
    --water-resname "$WATER_RESNAME" \
    --dist-cutoff "$HBOND_DIST_CUTOFF" \
    --angle-cutoff "$HBOND_ANGLE_CUTOFF" \
    --occ-cutoff "$WATER_CUTOFF_OCC" \
    --internal-residues $INTERNAL_RESIDS \
    --internal-water-ids $INTERNAL_WATER_IDS \
    --outdir "$RESULTS_DIR"

# ============================================================
# Step 10: Combined plots
# ============================================================
echo "=== Step 10: Combined plots ==="
python3 "$SCRIPTS_DIR/combined_plots.py" \
    --outdir "$RESULTS_DIR"

# ============================================================
# Step 11: Generate final Markdown report
# ============================================================
echo "=== Step 11: Generate final Markdown report ==="
python3 "$SCRIPTS_DIR/generate_report.py" \
    --outdir "$RESULTS_DIR"

# ============================================================
echo ""
echo "=============================================="
echo " Pipeline complete!"
echo " Report  : $RESULTS_DIR/analysis_report.md"
echo " Figures : $RESULTS_DIR/*.png"
echo " Tables  : $RESULTS_DIR/*.csv"
echo "=============================================="
