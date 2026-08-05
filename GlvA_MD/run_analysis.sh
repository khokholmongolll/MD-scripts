#!/bin/bash
set -e
# ============================================================
# Master Analysis Pipeline for MD Simulation (GlvA)
# Usage: ./path/to/run_analysis.sh [--force]
#
# Run from any directory — the script detects its own location
# to find the scripts/ folder.  All input/output paths are
# resolved relative to the current working directory unless
# given as absolute paths.
# ============================================================

# --------------- Locate this script ---------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$SCRIPT_DIR/scripts"

# --------------- Input files ---------------
# These are resolved relative to the CWD.  Use absolute paths if
# the files live outside the directory where you run the script.
TPR="step5_1.tpr"
GRO="step5_1.gro"
XTC_RAW="step5_1.xtc"
NDX_BASE="index_new.ndx"

# --------------- Directories ---------------
TMP_DIR="tmp"
RESULTS_DIR="results"

# --------------- Trajectory processing ---------------
# Point TRAJ_ANALYSIS to a pre-processed .xtc (absolute, or
# relative to CWD).  If you re-enable Steps 1-2 below, set
# TRAJ_ANALYSIS="$TRAJ_FIT" to use the freshly fitted trajectory.
TRAJ_PBC="$TMP_DIR/traj_pbc.xtc"
TRAJ_FIT="$TMP_DIR/4_traj_fit.xtc"
TRAJ_ANALYSIS="$TRAJ_FIT"

# --------------- Output files ---------------
NDX_UPDATED="$TMP_DIR/index_updated.ndx"
NDX_CMD_FILE="$TMP_DIR/ndx_commands.txt"

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

# --------------- Protein residues for MDAnalysis H-bond ---------------
PROTEIN_RESIDS="93 109 110 170 263 171 86 262 147 445 261 318 108"
LIGANDS_PHO="450"
LIGANDS_SUGARS="447 448 449"

# --------------- Water analysis ---------------
WATER_RESNAME="WAT"
WATER_CUTOFF_OCC=10.0
WATER_PROTEIN_RESIDS="93 109 110 170 263 108"

# Internal water analysis (specific water residue IDs)
INTERNAL_RESIDS="15 86 108 109 110 111 118 144 145 146 147 148 149 150 151 152 167 168 169 170 171 172 173 174 175 199 200 201 261 262 263 266 318 321 338 340 357"
INTERNAL_WATER_IDS="518 519 528 530 531 535 536 539 541 542 547 550 552 556 558 564 566 575"
INTERNAL_WATER_OCC_CUTOFF=5.0

# --------------- Bond distance atom pairs (0-based) ---------------
BOND_PAIRS=("2647 7039" "4108 7043")
BOND_LABELS=("ASH_170_OD2--4GA_448_O4" "TYD_263_OH--0GA_449_H2")

# --------------- RMSD skip groups ---------------
# Leave empty to compute all. Example: "System Water SOL non-Protein Water_and_ions K+ Cl-"
SKIP_RMSD_GROUPS="System Water SOL non-Protein Water_and_ions K+ Cl-"

# ============================================================
# Parse --force / -f flag
# ============================================================
FORCE=0
if [ "$1" = "-f" ] || [ "$1" = "--force" ]; then
    FORCE=1
    echo "[INFO] --force mode: all steps will execute regardless of cache."
fi

# ============================================================
# Caching infrastructure
# ============================================================

SKIPPED_STEPS=()
RAN_STEPS=()

# Compute a single md5 fingerprint from a list of dependencies.
# Each dependency is either a file path (gets md5summed) or a
# "name=value" string (included verbatim).
step_hash() {
    {
        for item in "$@"; do
            if [ -f "$item" ]; then
                printf "F:%s:%s\n" "$(basename "$item")" "$(md5sum "$item" | cut -d' ' -f1)"
            else
                printf "V:%s\n" "$item"
            fi
        done
    } | md5sum | cut -d' ' -f1
}

# Check whether a step can be skipped.
# Arguments: step_name <deps...>  --  <output_files...>
# Returns 0 if the step should be SKIPPED, 1 if it should RUN.
skip_step() {
    local step_name="$1"; shift

    # Separate dependencies from output files
    local deps=()
    local outputs=()
    local in_outputs=false
    for arg in "$@"; do
        if [ "$arg" = "--" ]; then
            in_outputs=true
        elif $in_outputs; then
            outputs+=("$arg")
        else
            deps+=("$arg")
        fi
    done

    # --force overrides caching → always run
    if [ "$FORCE" -eq 1 ]; then
        RAN_STEPS+=("${step_name}")
        return 1
    fi

    # All output files must exist
    for f in "${outputs[@]}"; do
        if [ ! -f "$f" ]; then
            echo "  [CACHE] Step ${step_name}: output '$f' missing, will run."
            RAN_STEPS+=("${step_name}")
            return 1
        fi
    done

    # Compute current fingerprint and compare with saved
    local state_file="$TMP_DIR/.step_${step_name}.state"
    local current_hash
    current_hash=$(step_hash "${deps[@]}")

    if [ -f "$state_file" ]; then
        local saved_hash
        saved_hash=$(cat "$state_file")
        if [ "$current_hash" = "$saved_hash" ]; then
            echo "  [CACHE] Step ${step_name}: inputs unchanged, skipping."
            SKIPPED_STEPS+=("${step_name}")
            return 0
        fi
    fi

    RAN_STEPS+=("${step_name}")
    return 1
}

# Write a step's state fingerprint after it completes.
save_step() {
    local step_name="$1"; shift
    local state_file="$TMP_DIR/.step_${step_name}.state"
    step_hash "$@" > "$state_file"
}

# ============================================================
# Prepare static dependency files
# ============================================================
mkdir -p "$TMP_DIR" "$RESULTS_DIR"

# ndx commands — always written so they can be hashed
cat > "$NDX_CMD_FILE" <<'EOF'
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

# Serialize arrays for hashing
BOND_DEPS_STR=""
for pair in "${BOND_PAIRS[@]}"; do BOND_DEPS_STR="${BOND_DEPS_STR}${pair}|"; done
BOND_DEPS_STR="${BOND_DEPS_STR}LABELS="
for lbl in "${BOND_LABELS[@]}"; do BOND_DEPS_STR="${BOND_DEPS_STR}${lbl}|"; done

# ============================================================
# STEP 1  —  PBC centering (commented by default)
# ============================================================
: <<'SKIP_STEP1'
echo "=== Step 1: PBC centering ==="
if skip_step 1 \
    "$TPR" "$XTC_RAW" "$NDX_BASE" \
    "CENTER_GROUP=$CENTER_GROUP" "OUTPUT_GROUP=$OUTPUT_GROUP" \
    -- "$TRAJ_PBC"
then
    :
else
    echo "$CENTER_GROUP $OUTPUT_GROUP" | gmx trjconv \
        -s "$TPR" -f "$XTC_RAW" -n "$NDX_BASE" \
        -pbc mol -center -ur compact \
        -o "$TRAJ_PBC"
    save_step 1 \
        "$TPR" "$XTC_RAW" "$NDX_BASE" \
        "CENTER_GROUP=$CENTER_GROUP" "OUTPUT_GROUP=$OUTPUT_GROUP"
fi
SKIP_STEP1

# ============================================================
# STEP 2  —  Fit rotation + translation (commented by default)
# ============================================================
: <<'SKIP_STEP2'
echo "=== Step 2: Fit rotation + translation ==="
if skip_step 2 \
    "$TPR" "$TRAJ_PBC" "$NDX_BASE" \
    "FIT_GROUP=$FIT_GROUP" "OUTPUT_GROUP=$OUTPUT_GROUP" \
    -- "$TRAJ_FIT"
then
    :
else
    echo "$FIT_GROUP $OUTPUT_GROUP" | gmx trjconv \
        -s "$TPR" -f "$TRAJ_PBC" -n "$NDX_BASE" \
        -fit rot+trans \
        -o "$TRAJ_FIT"
    save_step 2 \
        "$TPR" "$TRAJ_PBC" "$NDX_BASE" \
        "FIT_GROUP=$FIT_GROUP" "OUTPUT_GROUP=$OUTPUT_GROUP"
fi
SKIP_STEP2

# ============================================================
# STEP 3  —  Build custom index
# ============================================================
echo "=== Step 3: Building custom index ==="
if skip_step 3 \
    "$GRO" "$NDX_BASE" "$NDX_CMD_FILE" \
    -- "$NDX_UPDATED"
then
    :
else
    gmx make_ndx -f "$GRO" -n "$NDX_BASE" -o "$NDX_UPDATED" \
        < "$NDX_CMD_FILE" \
        > "$TMP_DIR/make_ndx.log" 2>&1
    echo "  Created $NDX_UPDATED"
    save_step 3 \
        "$GRO" "$NDX_BASE" "$NDX_CMD_FILE"
fi

# ============================================================
# STEP 4  —  GROMACS H-bond analysis (matrix)
# ============================================================
echo "=== Step 4: GROMACS H-bond analysis ==="
HB_OUT="$TMP_DIR/hbmat.xpm"
if skip_step 4 \
    "$TRAJ_ANALYSIS" "$TPR" "$NDX_UPDATED" \
    "GROUPS=$HBOND_GROUP_A,$HBOND_GROUP_B" \
    -- "$HB_OUT"
then
    :
else
    echo "$HBOND_GROUP_A $HBOND_GROUP_B" | gmx hbond \
        -f "$TRAJ_ANALYSIS" \
        -s "$TPR" \
        -n "$NDX_UPDATED" \
        -num "$TMP_DIR/hbnum.xvg" \
        -dist "$TMP_DIR/hbdist.xvg" \
        -hbm "$HB_OUT" \
        -hbn "$TMP_DIR/hbond.ndx" \
        -nomerge
    save_step 4 \
        "$TRAJ_ANALYSIS" "$TPR" "$NDX_UPDATED" \
        "GROUPS=$HBOND_GROUP_A,$HBOND_GROUP_B"
fi

# ============================================================
# STEP 5  —  H-bond matrix analysis + heatmap
# ============================================================
echo "=== Step 5: H-bond matrix analysis ==="
if skip_step 5 \
    "$GRO" "$TMP_DIR/hbond.ndx" "$HB_OUT" \
    "OCC_CUTOFF=$OCCUPANCY_CUTOFF" \
    -- "$RESULTS_DIR/hbond_matrix.png" "$RESULTS_DIR/hbond_matrix.csv"
then
    :
else
    python3 "$SCRIPTS_DIR/hbond_matrix.py" \
        --gro "$GRO" \
        --ndx "$TMP_DIR/hbond.ndx" \
        --xpm "$HB_OUT" \
        --cutoff "$OCCUPANCY_CUTOFF" \
        --outdir "$RESULTS_DIR"
    save_step 5 \
        "$GRO" "$TMP_DIR/hbond.ndx" "$HB_OUT" \
        "OCC_CUTOFF=$OCCUPANCY_CUTOFF"
fi

# ============================================================
# STEP 6  —  RMSD analysis for all index groups
# ============================================================
echo "=== Step 6: RMSD analysis ==="
if skip_step 6 \
    "$TPR" "$TRAJ_ANALYSIS" "$NDX_UPDATED" \
    "SKIP_RMSD=$SKIP_RMSD_GROUPS" \
    -- "$RESULTS_DIR/rmsd_overlay.png" "$RESULTS_DIR/rmsd_subplots.png" "$RESULTS_DIR/rmsd_data.csv"
then
    :
else
    RMSD_SKIP_ARGS=""
    [ -n "$SKIP_RMSD_GROUPS" ] && RMSD_SKIP_ARGS="--skip $SKIP_RMSD_GROUPS"
    python3 "$SCRIPTS_DIR/rmsd_analysis.py" \
        --tpr "$TPR" \
        --traj "$TRAJ_ANALYSIS" \
        --ndx "$NDX_UPDATED" \
        $RMSD_SKIP_ARGS \
        --outdir "$RESULTS_DIR"
    save_step 6 \
        "$TPR" "$TRAJ_ANALYSIS" "$NDX_UPDATED" \
        "SKIP_RMSD=$SKIP_RMSD_GROUPS"
fi

# ============================================================
# STEP 7  —  Bond distance analysis
# ============================================================
echo "=== Step 7: Bond distance analysis ==="
if skip_step 7 \
    "$TPR" "$TRAJ_ANALYSIS" \
    "BONDS=$BOND_DEPS_STR" \
    -- "$RESULTS_DIR/bond_distances.png" "$RESULTS_DIR/bond_distances.csv"
then
    :
else
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
    save_step 7 \
        "$TPR" "$TRAJ_ANALYSIS" \
        "BONDS=$BOND_DEPS_STR"
fi

# ============================================================
# STEP 8  —  Protein-ligand H-bond analysis
# ============================================================
echo "=== Step 8: Protein-ligand H-bond analysis ==="
if skip_step 8 \
    "$TPR" "$TRAJ_ANALYSIS" "$GRO" \
    "PROT_RESIDS=$PROTEIN_RESIDS" \
    "LIG_PHO=$LIGANDS_PHO" \
    "LIG_SUGARS=$LIGANDS_SUGARS" \
    "DIST_CUT=$HBOND_DIST_CUTOFF" \
    "ANGLE_CUT=$HBOND_ANGLE_CUTOFF" \
    "MIN_OCC=$MIN_OCCUPANCY_BONDS" \
    -- "$RESULTS_DIR/hbond_PHO.csv" "$RESULTS_DIR/hbond_4GA_0GA_ROH.csv"
then
    :
else
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
    save_step 8 \
        "$TPR" "$TRAJ_ANALYSIS" "$GRO" \
        "PROT_RESIDS=$PROTEIN_RESIDS" \
        "LIG_PHO=$LIGANDS_PHO" \
        "LIG_SUGARS=$LIGANDS_SUGARS" \
        "DIST_CUT=$HBOND_DIST_CUTOFF" \
        "ANGLE_CUT=$HBOND_ANGLE_CUTOFF" \
        "MIN_OCC=$MIN_OCCUPANCY_BONDS"
fi

# ============================================================
# STEP 9  —  Water-mediated H-bond analysis
# ============================================================
echo "=== Step 9: Water-mediated H-bond analysis ==="
if skip_step 9 \
    "$TPR" "$TRAJ_ANALYSIS" \
    "W_RESIDS=$WATER_PROTEIN_RESIDS" \
    "W_NAME=$WATER_RESNAME" \
    "W_DIST=$HBOND_DIST_CUTOFF" \
    "W_ANGLE=$HBOND_ANGLE_CUTOFF" \
    "W_OCC=$WATER_CUTOFF_OCC" \
    "I_RESIDS=$INTERNAL_RESIDS" \
    "I_WATERS=$INTERNAL_WATER_IDS" \
    "I_OCC=$INTERNAL_WATER_OCC_CUTOFF" \
    -- "$RESULTS_DIR/hbond_water.csv"
then
    :
else
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
        --internal-occ-cutoff "$INTERNAL_WATER_OCC_CUTOFF" \
        --outdir "$RESULTS_DIR"
    save_step 9 \
        "$TPR" "$TRAJ_ANALYSIS" \
        "W_RESIDS=$WATER_PROTEIN_RESIDS" \
        "W_NAME=$WATER_RESNAME" \
        "W_DIST=$HBOND_DIST_CUTOFF" \
        "W_ANGLE=$HBOND_ANGLE_CUTOFF" \
        "W_OCC=$WATER_CUTOFF_OCC" \
        "I_RESIDS=$INTERNAL_RESIDS" \
        "I_WATERS=$INTERNAL_WATER_IDS" \
        "I_OCC=$INTERNAL_WATER_OCC_CUTOFF"
fi

# ============================================================
# STEP 10  —  Combined plots
# ============================================================
echo "=== Step 10: Combined plots ==="
if skip_step 10 \
    "$RESULTS_DIR/rmsd_data.csv" "$RESULTS_DIR/bond_distances.csv" \
    -- "$RESULTS_DIR/combined_plots.png"
then
    :
else
    python3 "$SCRIPTS_DIR/combined_plots.py" \
        --outdir "$RESULTS_DIR"
    save_step 10 \
        "$RESULTS_DIR/rmsd_data.csv" "$RESULTS_DIR/bond_distances.csv"
fi

# ============================================================
# STEP 12  —  Metadynamics:  PLUMED sum_hills
# STEP 13  —  Metadynamics:  Python analysis + plots
#
# Runs only when HILLS and COLVAR exist in CWD.
# Skips gracefully otherwise (marker file for the report).
# ============================================================

META_FES2D="$TMP_DIR/fes2d.dat"
META_FES_PREFIX="$TMP_DIR/fes_"
META_BIN=100
META_STRIDE=1000
META_FES_PATTERN="${META_FES_PREFIX}[0-9]*"

if [ -f "HILLS" ] && [ -f "COLVAR" ]; then

    # --- 12: PLUMED ---
    echo "=== Step 12: Metadynamics — PLUMED sum_hills ==="
    if skip_step 12 \
        "HILLS" \
        "META_MIN=-pi,-pi" "META_MAX=pi,pi" \
        "META_BIN=$META_BIN" "META_STRIDE=$META_STRIDE" \
        -- "$META_FES2D"
    then
        :
    else
        plumed sum_hills --hills HILLS \
            --min -pi,-pi --max pi,pi --bin "$META_BIN,$META_BIN" \
            --outfile "$META_FES2D"
        plumed sum_hills --hills HILLS \
            --min -pi,-pi --max pi,pi --bin "$META_BIN,$META_BIN" \
            --stride "$META_STRIDE" --outfile "$META_FES_PREFIX"
        save_step 12 \
            "HILLS" \
            "META_MIN=-pi,-pi" "META_MAX=pi,pi" \
            "META_BIN=$META_BIN" "META_STRIDE=$META_STRIDE"
    fi

    # --- 13: Python analysis ---
    echo "=== Step 13: Metadynamics — Python analysis ==="
    META_CONV_LIST=$(find "$TMP_DIR" -maxdepth 1 -name 'fes_[0-9]*' \
                     -type f 2>/dev/null \
                     | sort | xargs echo 2>/dev/null || echo "")

    if skip_step 13 \
        "COLVAR" "$META_FES2D" \
        "CONV_FILES=$META_CONV_LIST" \
        -- "$RESULTS_DIR/meta_colvar.png" "$RESULTS_DIR/meta_fes2d.png"
    then
        :
    else
        python3 "$SCRIPTS_DIR/meta_analysis.py" \
            --colvar "COLVAR" \
            --fes2d "$META_FES2D" \
            --fes-pattern "$META_FES_PATTERN" \
            --outdir "$RESULTS_DIR"
        save_step 13 \
            "COLVAR" "$META_FES2D" \
            "CONV_FILES=$META_CONV_LIST"
    fi

    rm -f "$RESULTS_DIR/.meta_skipped"

else
    echo "=== Step 12-13: Metadynamics [SKIPPED] ==="
    echo "  HILLS and/or COLVAR not found in working directory."
    echo "  Place both files here to enable the metadynamics stage."
    touch "$RESULTS_DIR/.meta_skipped"
fi

# ============================================================
# STEP 11  —  Generate final Markdown report
# ============================================================
echo "=== Step 11: Generate final Markdown report ==="
# Depend on a snapshot of all results/ content
RESULTS_SNAPSHOT=$(find "$RESULTS_DIR" -type f \( -name "*.csv" -o -name "*.png" \) -exec md5sum {} \; 2>/dev/null | sort | md5sum | cut -d' ' -f1)

if skip_step 11 \
    "RESULTS_SNAPSHOT=$RESULTS_SNAPSHOT" \
    -- "$RESULTS_DIR/analysis_report.md"
then
    :
else
    python3 "$SCRIPTS_DIR/generate_report.py" \
        --outdir "$RESULTS_DIR"
    save_step 11 \
        "RESULTS_SNAPSHOT=$RESULTS_SNAPSHOT"
fi

# ============================================================
echo ""
echo "=============================================="
echo " Pipeline complete!"
echo " Skipped : ${#SKIPPED_STEPS[@]} steps  (${SKIPPED_STEPS[*]})"
echo " Run     : ${#RAN_STEPS[@]} steps  (${RAN_STEPS[*]})"
echo " Report  : $RESULTS_DIR/analysis_report.md"
echo " Figures : $RESULTS_DIR/*.png"
echo " Tables  : $RESULTS_DIR/*.csv"
echo "=============================================="
