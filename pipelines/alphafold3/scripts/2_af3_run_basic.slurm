#!/usr/bin/env bash
#SBATCH --job-name=af3_run
#SBATCH --partition=gpu_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --constraint=Milan|SapphireRapids
#SBATCH --mem=30G
#SBATCH --time=02:00:00

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

INPUT_LIST="${INPUT_LIST:-input.lst}"
JSON_DIR="${JSON_DIR:-jsons}"
OUT_DIR="${OUT_DIR:-outputs/af3_predictions/default_block}"
AF_OFFSET="${AF_OFFSET:-0}"

AF3_IMAGE="${AF3_IMAGE:-/apps/singularity-images/alphafold-3.0.0.sif}"
AF3_MODEL_DIR="${AF3_MODEL_DIR:-/scratch/mm07305/iob/af3_parameter}"
AF3_DB_DIR="${AF3_DB_DIR:-/db/AlphaFold3/20241114}"

[[ -f "$INPUT_LIST" ]] || { echo "ERROR: input list not found: $INPUT_LIST" >&2; exit 1; }
[[ -d "$JSON_DIR" ]] || { echo "ERROR: JSON directory not found: $JSON_DIR" >&2; exit 1; }
[[ -f "$AF3_IMAGE" ]] || { echo "ERROR: AF3 image not found: $AF3_IMAGE" >&2; exit 1; }
[[ -d "$AF3_MODEL_DIR" ]] || { echo "ERROR: AF3 model dir not found: $AF3_MODEL_DIR" >&2; exit 1; }
[[ -d "$AF3_DB_DIR" ]] || { echo "ERROR: AF3 DB dir not found: $AF3_DB_DIR" >&2; exit 1; }

task_id="${SLURM_ARRAY_TASK_ID:-1}"
line_no=$(( AF_OFFSET + task_id ))
json_file=$(awk -v n="$line_no" 'NR==n {print; exit}' "$INPUT_LIST")

[[ -n "${json_file:-}" ]] || { echo "ERROR: no JSON found at line ${line_no} in ${INPUT_LIST}" >&2; exit 1; }
[[ -f "${JSON_DIR}/${json_file}" ]] || { echo "ERROR: JSON file missing: ${JSON_DIR}/${json_file}" >&2; exit 1; }

mkdir -p "$OUT_DIR"

meta_file="${OUT_DIR}/task_${SLURM_JOB_ID}_${task_id}.meta.txt"
{
    echo "job_id=${SLURM_JOB_ID}"
    echo "task_id=${task_id}"
    echo "input_line=${line_no}"
    echo "json_file=${json_file}"
    echo "json_dir=${JSON_DIR}"
    echo "output_dir=${OUT_DIR}"
    echo "model_dir=${AF3_MODEL_DIR}"
    echo "db_dir=${AF3_DB_DIR}"
} > "$meta_file"

singularity exec \
  --nv \
  --bind "$PWD/${JSON_DIR}:/root/af_input" \
  --bind "$PWD/${OUT_DIR}:/root/af_output" \
  --bind "${AF3_MODEL_DIR}:/root/models" \
  --bind "${AF3_DB_DIR}:/root/public_databases" \
  "$AF3_IMAGE" \
  python /app/alphafold/run_alphafold.py \
    --json_path="/root/af_input/${json_file}" \
    --model_dir=/root/models \
    --db_dir=/root/public_databases \
    --jackhmmer_n_cpu "${SLURM_CPUS_PER_TASK}" \
    --output_dir=/root/af_output
