#!/bin/bash
#SBATCH --job-name=foldx_repair
#SBATCH --partition=bahl_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_foldx_repair_%j.out
#SBATCH --error=slurm_foldx_repair_%j.err

set -euo pipefail

module purge
export PATH="${HOME}/bin:${PATH}"

# ---------------- User-tunable settings ----------------
ION="${ION:-0.05}"
TEMP="${TEMP:-298}"
REPAIR_PH="${REPAIR_PH:-7.0}"
FOLDX_BIN="${FOLDX_BIN:-foldx}"
# ------------------------------------------------------

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch $0 <input.pdb>"
  echo "Optional env vars: ION=0.05 TEMP=298 REPAIR_PH=7.0 FOLDX_BIN=foldx"
  exit 1
fi

cd "${SLURM_SUBMIT_DIR:-$PWD}"

INPUT_PDB="$1"
[[ -f "$INPUT_PDB" ]] || { echo "ERROR: Input PDB not found: $INPUT_PDB"; exit 1; }

INPUT_ABS="$(python3 - <<'PY' "$INPUT_PDB"
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)"
INPUT_BASE="$(basename "$INPUT_PDB")"
INPUT_STEM="${INPUT_BASE%.pdb}"

OUTROOT="outputs/foldx_repair"
RUNROOT="runs/foldx_repair"
mkdir -p "$OUTROOT" "$RUNROOT"

WORK_PARENT="${SLURM_TMPDIR:-$RUNROOT}"
WORKDIR="${WORK_PARENT}/repair_${SLURM_JOB_ID:-manual}"
mkdir -p "$WORKDIR"
cp "$INPUT_ABS" "$WORKDIR/$INPUT_BASE"
cd "$WORKDIR"

echo "=========================================================="
echo " FoldX RepairPDB"
echo " Input:    $INPUT_BASE"
echo " pH:       $REPAIR_PH"
echo " Ion:      $ION"
echo " Temp:     $TEMP K"
echo " Workdir:  $WORKDIR"
echo "=========================================================="

"$FOLDX_BIN" --command RepairPDB \
  --pdb "$INPUT_BASE" \
  --pH "$REPAIR_PH" \
  --ionStrength "$ION" \
  --temperature "$TEMP"

REPAIRED_PDB="${INPUT_STEM}_Repair.pdb"
[[ -f "$REPAIRED_PDB" ]] || {
  echo "ERROR: Expected repaired PDB not found: $REPAIRED_PDB"
  ls -lh
  exit 1
}

DEST_DIR="${SLURM_SUBMIT_DIR:-$PWD}/$OUTROOT"
mkdir -p "$DEST_DIR"
cp -f "$REPAIRED_PDB" "$DEST_DIR/"
cp -f ./*.fxout "$DEST_DIR/" 2>/dev/null || true

printf 'input_pdb\trepaired_pdb\trepair_ph\tion_strength\ttemperature_K\n' > "$DEST_DIR/${INPUT_STEM}_repair_metadata.tsv"
printf '%s\t%s\t%s\t%s\t%s\n' \
  "$INPUT_BASE" "$REPAIRED_PDB" "$REPAIR_PH" "$ION" "$TEMP" >> "$DEST_DIR/${INPUT_STEM}_repair_metadata.tsv"

echo "[Done] Repaired model copied to: $OUTROOT/$REPAIRED_PDB"
