#!/bin/bash
#SBATCH --job-name=foldx_buildmodel_ph
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00
#SBATCH --array=1-1
#SBATCH --output=slurm_foldx_buildmodel_%A_%a.out
#SBATCH --error=slurm_foldx_buildmodel_%A_%a.err

set -euo pipefail
shopt -s nullglob

module purge
export PATH="${HOME}/bin:${PATH}"

# ---------------- User-tunable settings ----------------
MUTDIR="${MUTDIR:-foldx_mutfiles}"
NRUNS="${NRUNS:-3}"
ION="${ION:-0.05}"
TEMP="${TEMP:-298}"
PH_LIST_CSV="${PH_LIST_CSV:-4.5,5.0,5.3,5.5,5.7,6.0,6.5,7.0}"
FOLDX_BIN="${FOLDX_BIN:-foldx}"
# ------------------------------------------------------

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch --array=1-N $0 <repaired_input.pdb>"
  echo "Optional env vars: MUTDIR=foldx_mutfiles NRUNS=3 ION=0.05 TEMP=298 PH_LIST_CSV=4.5,5.0,..."
  exit 1
fi

cd "${SLURM_SUBMIT_DIR:-$PWD}"

REPAIRED_PDB="$1"
[[ -f "$REPAIRED_PDB" ]] || { echo "ERROR: Repaired PDB not found: $REPAIRED_PDB"; exit 1; }

REPAIRED_ABS="$(python3 - <<'PY' "$REPAIRED_PDB"
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)"
REPAIRED_BASE="$(basename "$REPAIRED_PDB")"
REPAIRED_STEM="${REPAIRED_BASE%.pdb}"
TEMP_TAG="T${TEMP}"

mapfile -t MUTFILES < <(
  find "$MUTDIR" -maxdepth 1 -type f \( -name '*.txt' -o -name '*.mut' \) \
    | sort -V \
    | grep -Ev 'individual_list_all|manifest'
)

MUT_COUNT="${#MUTFILES[@]}"
(( MUT_COUNT > 0 )) || { echo "ERROR: No mutation files found in $MUTDIR"; exit 1; }

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
(( TASK_ID >= 1 && TASK_ID <= MUT_COUNT )) || {
  echo "ERROR: SLURM_ARRAY_TASK_ID=$TASK_ID is outside 1..$MUT_COUNT"
  exit 1
}

MUTFILE_ABS="${MUTFILES[$((TASK_ID - 1))]}"
MUTFILE_BASE="$(basename "$MUTFILE_ABS")"
VARIANT_TAG="${MUTFILE_BASE%.*}"

OUTROOT="outputs/foldx_buildmodel/per_variant/${VARIANT_TAG}/${TEMP_TAG}"
FXOUT_DIR="$OUTROOT/fxout"
RUNROOT="runs/foldx_buildmodel/${VARIANT_TAG}/${TEMP_TAG}"
mkdir -p "$FXOUT_DIR" "$RUNROOT"

WORK_PARENT="${SLURM_TMPDIR:-$RUNROOT}"
WORKDIR="${WORK_PARENT}/job_${SLURM_ARRAY_JOB_ID:-manual}_${TASK_ID}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

cp "$REPAIRED_ABS" ./ref.pdb
tr -d '\r' < "$MUTFILE_ABS" > ./individual_list.txt
[[ -s individual_list.txt ]] || { echo "ERROR: Empty mutation file for $VARIANT_TAG"; exit 1; }

cp -f individual_list.txt "$OUTROOT/"
cp -f ref.pdb "$OUTROOT/${REPAIRED_BASE}"

IFS=',' read -r -a PH_LIST <<< "$PH_LIST_CSV"

FAILLOG="$OUTROOT/failures.tsv"
printf 'variant\ttemp\tpH\treturn_code\tnote\n' > "$FAILLOG"

printf 'variant\tmutfile\trepaired_pdb\ttemp\tnruns\tion_strength\tph_list\n' > "$OUTROOT/run_metadata.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$VARIANT_TAG" "$MUTFILE_BASE" "$REPAIRED_BASE" "$TEMP" "$NRUNS" "$ION" "$PH_LIST_CSV" >> "$OUTROOT/run_metadata.tsv"

echo "=========================================================="
echo " FoldX BuildModel pH scan"
echo " Variant:       $VARIANT_TAG"
echo " Mutation file: $MUTFILE_BASE"
echo " Repaired PDB:  $REPAIRED_BASE"
echo " Temp:          $TEMP K"
echo " pHs:           ${PH_LIST[*]}"
echo " Outdir:        $OUTROOT"
echo "=========================================================="

for PH in "${PH_LIST[@]}"; do
  PH_TAG="pH${PH//./p}"
  OUTTAG="${VARIANT_TAG}_${TEMP_TAG}_${PH_TAG}_${REPAIRED_STEM}"

  echo "[FoldX] TEMP=${TEMP} pH=${PH} -> ${OUTTAG}"

  set +e
  "$FOLDX_BIN" --command BuildModel \
    --pdb ref.pdb \
    --mutant-file individual_list.txt \
    --numberOfRuns "$NRUNS" \
    --pH "$PH" \
    --ionStrength "$ION" \
    --temperature "$TEMP" \
    --out-pdb false \
    --output-dir "$WORKDIR" \
    --output-file "$OUTTAG"
  rc=$?
  set -e

  if (( rc != 0 )); then
    printf '%s\t%s\t%s\t%s\t%s\n' "$VARIANT_TAG" "$TEMP" "$PH" "$rc" 'foldx_failed' >> "$FAILLOG"
    echo "  [WARN] FoldX failed (rc=$rc) at pH=${PH}; continuing"
  fi

  for f in Dif_"${OUTTAG}"*.fxout Average_"${OUTTAG}"*.fxout Raw_"${OUTTAG}"*.fxout PdbList_"${OUTTAG}"*.fxout ; do
    [[ -f "$f" ]] && mv -f "$f" "$FXOUT_DIR/"
  done
done

echo "[Done] $VARIANT_TAG -> $OUTROOT"
echo "  - FoldX outputs: $FXOUT_DIR"
echo "  - Failure log:   $FAILLOG"
