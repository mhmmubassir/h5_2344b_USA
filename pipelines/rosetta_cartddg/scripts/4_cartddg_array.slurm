#!/bin/bash
#SBATCH --job-name=cartddg
#SBATCH --partition=bahl_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=100:00:00
#SBATCH --array=1-999
#SBATCH --output=slurm_cartddg_%A_%a.out
#SBATCH --error=slurm_cartddg_%A_%a.err

set -euo pipefail

ROS_MODULE="${ROS_MODULE:-Rosetta/2022.46.334-intel-2021b}"
ROSETTA_DB="${ROSETTA_DB:-/apps/eb/Rosetta/2022.46.334-intel-2021b/database}"
MUT_DIR="${MUT_DIR:-inputs/mutfiles}"
ITERATIONS="${ITERATIONS:-5}"
BBNBRS="${BBNBRS:-1}"

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch --array=1-N $0 <relaxed_input.pdb>"
    exit 1
fi

cd "$SLURM_SUBMIT_DIR" || {
    echo "ERROR: cannot cd to SLURM_SUBMIT_DIR=$SLURM_SUBMIT_DIR"
    exit 1
}

INPUT_PDB="$1"
[[ -f "$INPUT_PDB" ]] || {
    echo "ERROR: input PDB not found: $INPUT_PDB"
    exit 1
}

module load "$ROS_MODULE"

mapfile -t MUT_FILES < <(find "$MUT_DIR" -maxdepth 1 -type f -name '*.mut' | sort -V)
(( ${#MUT_FILES[@]} > 0 )) || {
    echo "ERROR: no .mut files found in $MUT_DIR"
    exit 1
}

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
if (( TASK_ID < 1 || TASK_ID > ${#MUT_FILES[@]} )); then
    echo "Array index $TASK_ID is outside the discovered mut-file range 1..${#MUT_FILES[@]}"
    echo "Tip: submit with a matching array range, for example --array=1-${#MUT_FILES[@]}"
    exit 1
fi

MUT_FILE="${MUT_FILES[$((TASK_ID - 1))]}"
MUT_TAG=$(basename "$MUT_FILE" .mut)
INPUT_BASE=$(basename "$INPUT_PDB")
INPUT_STEM="${INPUT_BASE%.pdb}"

RUN_DIR="runs/cartddg/${MUT_TAG}"
DDG_DIR="outputs/cartddg/ddg"
SCORE_DIR="outputs/cartddg/scores"
WT_DIR="outputs/cartddg/wt_pdbs"
MUT_PDB_DIR="outputs/cartddg/mut_pdbs"
LOG_DIR="logs/cartddg"

mkdir -p "$RUN_DIR" "$DDG_DIR" "$SCORE_DIR" "$WT_DIR" "$MUT_PDB_DIR" "$LOG_DIR"

cp "$INPUT_PDB" "$RUN_DIR/"
cp "$MUT_FILE" "$RUN_DIR/"

cd "$RUN_DIR"

LOCAL_PDB=$(basename "$INPUT_PDB")
LOCAL_MUT=$(basename "$MUT_FILE")

echo "============================================================"
echo "Cartesian ddG job"
echo "Input relaxed PDB : $INPUT_PDB"
echo "Mutation file     : $MUT_FILE"
echo "Mutation tag      : $MUT_TAG"
echo "Array index       : $TASK_ID / ${#MUT_FILES[@]}"
echo "Run directory     : $RUN_DIR"
echo "============================================================"

srun --mpi=pmix_v3 cartesian_ddg.mpi.linuxiccrelease \
  -database "$ROSETTA_DB" \
  -in:file:s "$LOCAL_PDB" \
  -ddg::mut_file "$LOCAL_MUT" \
  -ddg::cartesian \
  -ddg::dump_pdbs true \
  -ddg::legacy false \
  -ddg::iterations "$ITERATIONS" \
  -ddg::bbnbrs "$BBNBRS" \
  -ddg::min_cst true \
  -fa_max_dis 9.0 \
  -beta_nov16_cart \
  -ideal_sugars \
  -auto_detect_glycan_connections \
  -maintain_links \
  -include_sugars \
  -alternate_3_letter_codes pdb_sugar \
  -ignore_unrecognized_res false \
  -ignore_zero_occupancy false \
  -write_pdb_link_records \
  -ex1 -ex2

shopt -s nullglob

for wt_pdb in WT_*.pdb; do
    mv "$wt_pdb" "$SLURM_SUBMIT_DIR/$WT_DIR/WT_${MUT_TAG}_${wt_pdb#WT_}"
done

for mut_pdb in MUT_*.pdb; do
    mv "$mut_pdb" "$SLURM_SUBMIT_DIR/$MUT_PDB_DIR/MUT_${MUT_TAG}_${mut_pdb#MUT_}"
done

for ddgfile in *.ddg; do
    mv "$ddgfile" "$SLURM_SUBMIT_DIR/$DDG_DIR/${MUT_TAG}_$(basename "$ddgfile")"
done

for scfile in *.sc; do
    mv "$scfile" "$SLURM_SUBMIT_DIR/$SCORE_DIR/${MUT_TAG}_$(basename "$scfile")"
done

OUT_FILE="$SLURM_SUBMIT_DIR/slurm_cartddg_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
ERR_FILE="$SLURM_SUBMIT_DIR/slurm_cartddg_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
[[ -f "$OUT_FILE" ]] && cp "$OUT_FILE" "$LOG_DIR/${MUT_TAG}.out" && rm -f "$OUT_FILE"
[[ -f "$ERR_FILE" ]] && cp "$ERR_FILE" "$LOG_DIR/${MUT_TAG}.err" && rm -f "$ERR_FILE"

echo "Done. Results for $MUT_TAG are in:"
echo "  $DDG_DIR"
echo "  $SCORE_DIR"
echo "  $WT_DIR"
echo "  $MUT_PDB_DIR"
