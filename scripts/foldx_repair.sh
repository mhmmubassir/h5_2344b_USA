#!/bin/bash
#SBATCH --job-name=foldx_repair
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=02:00:00
#SBATCH --output=logs/foldx_repair_%j.out
#SBATCH --error=logs/foldx_repair_%j.err
#SBATCH --mail-user=mm07305@uga.edu
#SBATCH --mail-type=ALL

set -euo pipefail
module purge
export PATH="$HOME/bin:$PATH"

# ---- EDIT if needed ----
PDB_IN="MUT_EPI_ISL_18133029_4muts.pdb"
ION=0.05
TEMP=298
REPAIR_PH=7.0
# ------------------------

# IMPORTANT: logs/ must exist BEFORE sbatch is submitted (Slurm opens log files early)
# mkdir -p logs  # (keep for safety, but create it in submit dir before sbatch)
mkdir -p repaired

cd "$SLURM_SUBMIT_DIR"

if [[ ! -f "$PDB_IN" ]]; then
  echo "ERROR: PDB not found in submit dir: $PDB_IN"
  exit 1
fi

WORK="/scratch/mm07305/iob/foldx_jobs/repair_${SLURM_JOB_ID}"
mkdir -p "$WORK"
cp "$PDB_IN" "$WORK/"
cd "$WORK"

echo "[FoldX] RepairPDB on: $PDB_IN  (pH=${REPAIR_PH}, I=${ION}, T=${TEMP}K)"
foldx --command RepairPDB \
      --pdb "$PDB_IN" \
      --pH "$REPAIR_PH" \
      --ionStrength "$ION" \
      --temperature "$TEMP"

REPAIRED="${PDB_IN%.pdb}_Repair.pdb"
if [[ ! -f "$REPAIRED" ]]; then
  echo "ERROR: Expected repaired PDB not found: $REPAIRED"
  ls -lh
  exit 1
fi

cp -a "$REPAIRED" "$SLURM_SUBMIT_DIR/repaired/"
# keep the FoldX report too (handy for audit)
cp -a *.fxout "$SLURM_SUBMIT_DIR/repaired/" 2>/dev/null || true

echo "[Done] Repaired PDB written to: repaired/$REPAIRED"
