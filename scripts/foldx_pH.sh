#!/bin/bash
#SBATCH --job-name=foldx_pH_T298
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00
#SBATCH --array=1-2116
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -u
set -o pipefail

module purge
export PATH="$HOME/bin:$PATH"

# ---------------- EDIT IF NEEDED ----------------
PDB_REPAIRED="repaired/MUT_EPI_ISL_18133029_4muts_Repair.pdb"
MUTDIR="foldx_mutfiles"
NRUNS=3
ION=0.05
TEMP=298

# wetlab-like anchors included
PH_LIST=(4.5 5.0 5.3 5.5 5.7 6.0 6.5 7.0)
# ------------------------------------------------

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs

OUTROOT="$SLURM_SUBMIT_DIR/foldx_ph_out/per_variant"
mkdir -p "$OUTROOT"

PDB_ABS="$SLURM_SUBMIT_DIR/$PDB_REPAIRED"
[[ -f "$PDB_ABS" ]] || { echo "ERROR: Repaired PDB not found: $PDB_ABS"; exit 1; }

# Collect mutfiles in stable order
shopt -s nullglob
candidates=( "$SLURM_SUBMIT_DIR/${MUTDIR}"/*__individual_list_[0-9]*.txt \
             "$SLURM_SUBMIT_DIR/${MUTDIR}"/individual_list_[0-9]*.txt )
filtered=()
for f in "${candidates[@]}"; do
  [[ "$(basename "$f")" == "individual_list_all.txt" ]] && continue
  filtered+=( "$f" )
done
mapfile -t mutation_files < <(printf '%s\n' "${filtered[@]}" | sort)

mut_count=${#mutation_files[@]}
(( mut_count > 0 )) || { echo "ERROR: No per-variant mutation files in $MUTDIR"; exit 1; }
(( SLURM_ARRAY_TASK_ID <= mut_count )) || { echo "ERROR: SLURM_ARRAY_TASK_ID too large"; exit 1; }

mut_file_abs="${mutation_files[$SLURM_ARRAY_TASK_ID - 1]}"
mut_base="$(basename "$mut_file_abs")"
mut_prefix="${mut_base%.txt}"

[[ -f "$mut_file_abs" ]] || { echo "ERROR: mutfile missing: $mut_file_abs"; exit 1; }

VDIR="$OUTROOT/$mut_prefix"
mkdir -p "$VDIR"
cd "$VDIR"

# Keep inputs next to outputs
cp -f "$PDB_ABS" ./ref.pdb
tr -d '\r' < "$mut_file_abs" > ./individual_list.txt
[[ -s individual_list.txt ]] || { echo "ERROR: empty individual_list.txt for $mut_prefix"; exit 1; }

FAILLOG="failures.tsv"
echo -e "variant\tpH\treturn_code\tnote" > "$FAILLOG"

echo "=========================================================="
echo " Job:     ${SLURM_ARRAY_JOB_ID:-NA}"
echo " Task:    ${SLURM_ARRAY_TASK_ID} / ${mut_count}"
echo " Variant: $mut_prefix"
echo " Outdir:  $VDIR"
echo " TEMP:    $TEMP K"
echo " pHs:     ${PH_LIST[*]}"
echo "=========================================================="

mkdir -p fxout

for PH in "${PH_LIST[@]}"; do
  TAG="pH$(echo "$PH" | sed 's/\./p/g')"
  OUTTAG="${mut_prefix}_T${TEMP}_${TAG}_ref"

  echo "[FoldX] BuildModel TEMP=${TEMP} pH=${PH} -> ${OUTTAG}"

  # Continue even if FoldX fails at one pH
  set +e
  foldx --command BuildModel \
        --pdb "ref.pdb" \
        --mutant-file "individual_list.txt" \
        --numberOfRuns "$NRUNS" \
        --pH "$PH" \
        --ionStrength "$ION" \
        --temperature "$TEMP" \
        --out-pdb false \
        --output-dir "$VDIR" \
        --output-file "$OUTTAG"
  rc=$?
  set -e

  if (( rc != 0 )); then
    echo -e "${mut_prefix}\t${PH}\t${rc}\tfoldx_failed" >> "$FAILLOG"
    echo "  [WARN] FoldX failed (rc=$rc) at pH=${PH}; continuing..."
  fi

  # Move main outputs into fxout/ to keep folder clean
  for f in Dif_"${OUTTAG}"*.fxout Average_"${OUTTAG}"*.fxout Raw_"${OUTTAG}"*.fxout PdbList_"${OUTTAG}"*.fxout; do
    [[ -f "$f" ]] && mv -f "$f" fxout/
  done

  # Optional: FoldX sometimes creates molecules/
  # Keep it if you want; otherwise uncomment:
  # rm -rf molecules 2>/dev/null || true
done

echo "[Done] $mut_prefix -> $VDIR"
echo "  - outputs: $VDIR/fxout/"
echo "  - failures: $VDIR/$FAILLOG"
