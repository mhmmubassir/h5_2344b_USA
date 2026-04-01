#!/bin/bash
#SBATCH --job-name=cart_relax
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_cart_relax_%j.out
#SBATCH --error=slurm_cart_relax_%j.err

set -euo pipefail

# Optional overrides at submission time, for example:
#   sbatch --export=NSTRUCT=10,RELAX_PARTITION=bahl_p 1_cart_relax.slurm input.pdb
ROS_MODULE="${ROS_MODULE:-Rosetta/2022.46.334-intel-2021b}"
ROSETTA_DB="${ROSETTA_DB:-/apps/eb/Rosetta/2022.46.334-intel-2021b/database}"
RELAX_SCRIPT="${RELAX_SCRIPT:-protocols/cart_relax.script}"
NSTRUCT="${NSTRUCT:-20}"

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch $0 <input_glycoprotein.pdb>"
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

[[ -f "$RELAX_SCRIPT" ]] || {
    echo "ERROR: relax script not found: $RELAX_SCRIPT"
    exit 1
}

module load "$ROS_MODULE"

INPUT_BASE=$(basename "$INPUT_PDB")
INPUT_STEM="${INPUT_BASE%.pdb}"
RUN_DIR="runs/cart_relax/${INPUT_STEM}"
MODEL_DIR="outputs/cart_relax/models/${INPUT_STEM}"
SCORE_DIR="outputs/cart_relax/scores"
BEST_DIR="outputs/cart_relax/best_models"
LOG_DIR="logs/cart_relax"

mkdir -p "$RUN_DIR" "$MODEL_DIR" "$SCORE_DIR" "$BEST_DIR" "$LOG_DIR"

cp "$INPUT_PDB" "$RUN_DIR/"
cp "$RELAX_SCRIPT" "$RUN_DIR/cart_relax.script"

cd "$RUN_DIR"

echo "============================================================"
echo "Cartesian relax job"
echo "Input PDB      : $INPUT_BASE"
echo "Input stem     : $INPUT_STEM"
echo "nstruct        : $NSTRUCT"
echo "Run directory  : $RUN_DIR"
echo "Model output   : $MODEL_DIR"
echo "Score output   : $SCORE_DIR/${INPUT_STEM}_relax.sc"
echo "============================================================"

srun --mpi=pmix_v3 relax.mpi.linuxiccrelease \
  -database "$ROSETTA_DB" \
  -s "$INPUT_BASE" \
  -use_input_sc \
  -nstruct "$NSTRUCT" \
  -out:file:scorefile score.sc \
  -relax:cartesian \
  -relax:script cart_relax.script \
  -relax:min_type lbfgs_armijo_nonmonotone \
  -beta_nov16_cart \
  -fa_max_dis 9.0 \
  -constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
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
for pdb in *.pdb; do
    cp "$pdb" "$MODEL_DIR/"
done
cp score.sc "$SCORE_DIR/${INPUT_STEM}_relax.sc"

OUT_FILE="$SLURM_SUBMIT_DIR/slurm_cart_relax_${SLURM_JOB_ID}.out"
ERR_FILE="$SLURM_SUBMIT_DIR/slurm_cart_relax_${SLURM_JOB_ID}.err"
[[ -f "$OUT_FILE" ]] && cp "$OUT_FILE" "$LOG_DIR/${INPUT_STEM}.out" && rm -f "$OUT_FILE"
[[ -f "$ERR_FILE" ]] && cp "$ERR_FILE" "$LOG_DIR/${INPUT_STEM}.err" && rm -f "$ERR_FILE"

echo "Done. Relax models are in $MODEL_DIR"
echo "Score file is in $SCORE_DIR/${INPUT_STEM}_relax.sc"
echo "Next steps:"
echo "  python scripts/2_sort_relax_scores.py $SCORE_DIR/${INPUT_STEM}_relax.sc"
echo "  python scripts/3_select_best_relaxed_model.py $SCORE_DIR/${INPUT_STEM}_relax.sc $MODEL_DIR"
