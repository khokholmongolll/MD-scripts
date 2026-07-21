#!/bin/bash
# ============================================================
# GROMACS + PLUMED Metadynamics Launch Script (Non-SLURM)
# ============================================================
# Usage: ./run_metad.sh
# Assumes GROMACS and PLUMED are installed and sourced in environment.
# Edit the CONFIGURATION section below for your specific system.
# ============================================================

set -euo pipefail  # Exit on error, undefined var, pipe failure

# ============================================================
# CONFIGURATION - EDIT THESE FOR YOUR SYSTEM
# ============================================================
# Project name (used for output files)
PROJECT="glvA_metad"

# Input files (relative to this script's directory)
GRO_FILE="step5_1.gro"          # Starting coordinates
TOP_FILE="topol.top"            # Topology
MDP_FILE="md.mdp"               # MD parameters (nvt/npt/prod)
PLUMED_FILE="plumed.dat"        # PLUMED input (this directory)

# GROMACS executable (if not in PATH, set full path or 'gmx_mpi')
GMX="gmx"

# Simulation control
MAX_WARN=1                      # Max warnings for grompp
NT_THREADS=0                    # OpenMP threads (0 = auto)
NT_MPI=1                        # MPI ranks (1 for single GPU/CPU)
GPU_ID="0"                      # GPU ID (if using GPU)
PIN_OFFSET=0                    # Thread pinning offset

# Output control
LOG_DIR="logs"
mkdir -p "${LOG_DIR}"

# ============================================================
# DERIVED SETTINGS - USUALLY NO NEED TO CHANGE
# ============================================================
TPR_FILE="${PROJECT}.tpr"
TRR_FILE="${PROJECT}.trr"
XTC_FILE="${PROJECT}.xtc"
EDR_FILE="${PROJECT}.edr"
LOG_FILE="${LOG_DIR}/${PROJECT}.log"
HILLS_FILE="HILLS"              # PLUMED hills output
COLVAR_FILE="COLVAR"            # PLUMED colvar output

# ============================================================
# PRE-FLIGHT CHECKS
# ============================================================
echo "=== PRE-FLIGHT CHECKS ==="
for f in "${GRO_FILE}" "${TOP_FILE}" "${MDP_FILE}" "${PLUMED_FILE}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Required input file not found: ${f}" >&2
        exit 1
    fi
done

if ! command -v "${GMX}" &> /dev/null; then
    echo "ERROR: GROMACS executable '${GMX}' not found in PATH." >&2
    echo "       Source your GMXRC or set GMX variable." >&2
    exit 1
fi

# Check PLUMED is linked/working
if ! "${GMX}" mdrun -h 2>&1 | grep -q "plumed"; then
    echo "WARNING: PLUMED support not detected in '${GMX} mdrun'." >&2
    echo "         Ensure GROMACS was compiled with PLUMED." >&2
fi

# ============================================================
# GOMPP - Generate TPR
# ============================================================
echo "=== RUNNING GOMPP ==="
"${GMX}" grompp \
    -f "${MDP_FILE}" \
    -c "${GRO_FILE}" \
    -p "${TOP_FILE}" \
    -o "${TPR_FILE}" \
    -maxwarn "${MAX_WARN}" \
    2>&1 | tee "${LOG_DIR}/grompp.log"

if [[ ! -f "${TPR_FILE}" ]]; then
    echo "ERROR: grompp failed to produce ${TPR_FILE}" >&2
    exit 1
fi

# ============================================================
# MDRUN - Production Metadynamics
# ============================================================
echo "=== STARTING METADYNAMICS (MDRUN) ==="
echo "Outputs:"
echo "  TPR:      ${TPR_FILE}"
echo "  TRR:      ${TRR_FILE}"
echo "  XTC:      ${XTC_FILE}"
echo "  EDR:      ${EDR_FILE}"
echo "  LOG:      ${LOG_FILE}"
echo "  HILLS:    ${HILLS_FILE}"
echo "  COLVAR:   ${COLVAR_FILE}"
echo ""

# Build mdrun command
MDRUN_CMD=(
    "${GMX}" mdrun
    -s "${TPR_FILE}"
    -o "${TRR_FILE}"
    -x "${XTC_FILE}"
    -e "${EDR_FILE}"
    -g "${LOG_FILE}"
    -plumed "${PLUMED_FILE}"
    -ntomp "${NT_THREADS}"
    -ntmpi "${NT_MPI}"
    -pin on
    -pinoffset "${PIN_OFFSET}"
)

# Add GPU flags if GPU_ID is set and non-empty
if [[ -n "${GPU_ID}" ]]; then
    MDRUN_CMD+=(-gpu_id "${GPU_ID}")
fi

# Execute
echo "Command: ${MDRUN_CMD[*]}"
"${MDRUN_CMD[@]}" 2>&1 | tee -a "${LOG_FILE}"

# ============================================================
# POST-RUN CHECK
# ============================================================
echo "=== SIMULATION FINISHED ==="
if [[ -f "${HILLS_FILE}" && -f "${COLVAR_FILE}" ]]; then
    echo "SUCCESS: PLUMED output files found."
    echo "  Hills lines: $(wc -l < "${HILLS_FILE}")"
    echo "  Colvar lines: $(wc -l < "${COLVAR_FILE}")"
else
    echo "WARNING: PLUMED output files (HILLS, COLVAR) not found." >&2
    echo "         Check ${LOG_FILE} for errors." >&2
fi

# Quick energy check
if [[ -f "${EDR_FILE}" ]]; then
    echo ""
    echo "=== ENERGY SUMMARY (last frame) ==="
    echo "Potential  Kinetic  Total  Temp  Pressure" | column -t
    "${GMX}" energy -f "${EDR_FILE}" -o /dev/stdout <<< "Potential Kinetic Total Temperature Pressure" 2>/dev/null | tail -1 | column -t
fi

echo "Done."