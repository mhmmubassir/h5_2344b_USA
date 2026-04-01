#!/usr/bin/env bash
set -euo pipefail

INPUT_LIST="${INPUT_LIST:-input.lst}"
BLOCK_SIZE="${BLOCK_SIZE:-8}"
RUNNER="${RUNNER:-2_af3_run_basic.slurm}"
JSON_DIR="${JSON_DIR:-jsons}"
OUTPUT_ROOT="${OUTPUT_ROOT:-outputs/af3_predictions}"
LOG_ROOT="${LOG_ROOT:-logs/af3}"

[[ -f "$INPUT_LIST" ]] || { echo "ERROR: input list not found: $INPUT_LIST" >&2; exit 1; }
[[ -f "$RUNNER" ]] || { echo "ERROR: runner script not found: $RUNNER" >&2; exit 1; }
[[ -d "$JSON_DIR" ]] || { echo "ERROR: JSON directory not found: $JSON_DIR" >&2; exit 1; }

mkdir -p "$OUTPUT_ROOT" "$LOG_ROOT"

TOTAL=$(wc -l < "$INPUT_LIST")
[[ "$TOTAL" -ge 1 ]] || { echo "ERROR: input list is empty: $INPUT_LIST" >&2; exit 1; }

offset=0
block_num=1

while [[ "$offset" -lt "$TOTAL" ]]; do
    remain=$(( TOTAL - offset ))
    tasks=$(( remain < BLOCK_SIZE ? remain : BLOCK_SIZE ))
    range="1-${tasks}"

    block_id=$(printf "block%04d" "$block_num")
    out_dir="${OUTPUT_ROOT}/${block_id}"
    log_dir="${LOG_ROOT}/${block_id}"
    mkdir -p "$out_dir" "$log_dir"

    echo "Submitting ${block_id} for input lines $((offset + 1))-$((offset + tasks))"

    jobid=$(sbatch --parsable \
        --array="$range" \
        --export=ALL,AF_OFFSET="$offset",INPUT_LIST="$INPUT_LIST",JSON_DIR="$JSON_DIR",OUT_DIR="$out_dir" \
        --output="${log_dir}/af3_%A_%a.out" \
        --error="${log_dir}/af3_%A_%a.err" \
        "$RUNNER")

    echo "  queued job ${jobid}"

    while squeue -h -j "$jobid" | grep -q .; do
        sleep 120
    done

    offset=$(( offset + tasks ))
    block_num=$(( block_num + 1 ))
done

echo "All ${TOTAL} AF3 inputs submitted and completed."
