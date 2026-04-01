#!/usr/bin/env bash
#SBATCH --job-name=flexddg
#SBATCH --partition=bahl_p
#SBATCH --constraint=EDR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=96:00:00
#SBATCH --array=1-1
#SBATCH --output=slurm_flexddg_%A_%a.out
#SBATCH --error=slurm_flexddg_%A_%a.err

# -----------------------------------------------------------------------------
# Generic SLURM array launcher for Rosetta Flex ddG with glycans.
#
# Usage:
#   sbatch --array=1-N run_flexddg_array.slurm <pre_minimized_complex.pdb>
#
# Expected repo layout:
#   protocols/ddG-backrub.xml
#   inputs/resfiles/*.resfile
#   outputs/{ddg_db3,score_sc,struct_db3,pdb_snapshots}
#
# Optional environment-variable overrides:
#   ROSETTA_MODULE=Rosetta/2022.46.334-intel-2021b
#   ROSETTA_DB=/apps/eb/Rosetta/2022.46.334-intel-2021b/database
#   RESFILE_DIR=inputs/resfiles
#   PROTOCOL_XML=protocols/ddG-backrub.xml
#   OUTPUT_ROOT=outputs
#   RUN_ROOT=runs
# -----------------------------------------------------------------------------

set -euo pipefail
shopt -s nullglob

readonly ROSETTA_MODULE="${ROSETTA_MODULE:-Rosetta/2022.46.334-intel-2021b}"
readonly ROSETTA_DB="${ROSETTA_DB:-/apps/eb/Rosetta/2022.46.334-intel-2021b/database}"
readonly RESFILE_DIR="${RESFILE_DIR:-inputs/resfiles}"
readonly PROTOCOL_XML="${PROTOCOL_XML:-protocols/ddG-backrub.xml}"
readonly OUTPUT_ROOT="${OUTPUT_ROOT:-outputs}"
readonly RUN_ROOT="${RUN_ROOT:-runs}"
readonly SNAPSHOT_DIR="${OUTPUT_ROOT}/pdb_snapshots"
readonly SCORE_DIR="${OUTPUT_ROOT}/score_sc"
readonly DDG_DIR="${OUTPUT_ROOT}/ddg_db3"
readonly STRUCT_DIR="${OUTPUT_ROOT}/struct_db3"

usage() {
  echo "Usage: sbatch --array=1-N $0 <pre_minimized_complex.pdb>"
}

die() {
  echo "❌ $*" >&2
  exit 1
}

info() {
  echo "[$(date '+%F %T')] $*"
}

if [[ $# -ne 1 ]]; then
  usage
  exit 1
fi

INPUT_PDB=$1
[[ -f "$INPUT_PDB" ]] || die "Input PDB not found: $INPUT_PDB"
[[ -f "$PROTOCOL_XML" ]] || die "Protocol XML not found: $PROTOCOL_XML"
[[ -d "$RESFILE_DIR" ]] || die "Resfile directory not found: $RESFILE_DIR"

cd "${SLURM_SUBMIT_DIR:-$PWD}"

mkdir -p "$RUN_ROOT" "$SNAPSHOT_DIR" "$SCORE_DIR" "$DDG_DIR" "$STRUCT_DIR"

module load "$ROSETTA_MODULE"
command -v rosetta_scripts.mpi.linuxiccrelease >/dev/null 2>&1 \
  || die "rosetta_scripts.mpi.linuxiccrelease not found after loading $ROSETTA_MODULE"
[[ -d "$ROSETTA_DB" ]] || die "Rosetta database directory not found: $ROSETTA_DB"

mapfile -t RESFILES < <(find "$RESFILE_DIR" -maxdepth 1 -type f -name '*.resfile' | sort -V)
(( ${#RESFILES[@]} > 0 )) || die "No .resfile files found in $RESFILE_DIR"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
[[ "$TASK_ID" =~ ^[0-9]+$ ]] || die "Invalid SLURM_ARRAY_TASK_ID: $TASK_ID"
(( TASK_ID >= 1 && TASK_ID <= ${#RESFILES[@]} )) || \
  die "Array index $TASK_ID is outside 1-${#RESFILES[@]}"

RESFILE="${RESFILES[TASK_ID-1]}"
MUT_TAG="$(basename "$RESFILE" .resfile)"
RUN_DIR="$RUN_ROOT/$MUT_TAG"
INPUT_PDB_BASENAME="$(basename "$INPUT_PDB")"
RESFILE_BASENAME="$(basename "$RESFILE")"
PROTOCOL_BASENAME="$(basename "$PROTOCOL_XML")"

mkdir -p "$RUN_DIR"
cp -f "$INPUT_PDB" "$RUN_DIR/$INPUT_PDB_BASENAME"
cp -f "$RESFILE" "$RUN_DIR/$RESFILE_BASENAME"
cp -f "$PROTOCOL_XML" "$RUN_DIR/$PROTOCOL_BASENAME"

cat > "$RUN_DIR/run_metadata.txt" <<META
mutation_tag=$MUT_TAG
array_index=$TASK_ID
input_pdb=$INPUT_PDB_BASENAME
resfile=$RESFILE_BASENAME
protocol_xml=$PROTOCOL_BASENAME
rosetta_module=$ROSETTA_MODULE
rosetta_database=$ROSETTA_DB
META

cd "$RUN_DIR"

info "Starting Flex ddG for mutation: $MUT_TAG"
info "Array index: $TASK_ID / ${#RESFILES[@]}"

srun --mpi=pmix_v3 -n "${SLURM_NTASKS:-1}" \
  rosetta_scripts.mpi.linuxiccrelease \
  -database "$ROSETTA_DB" \
  -in:file:s "$INPUT_PDB_BASENAME" \
  -in:file:fullatom \
  -parser:protocol "$PROTOCOL_BASENAME" \
  -parser:script_vars \
      chainstomove=A \
      mutate_resfile_relpath="$RESFILE_BASENAME" \
      number_backrub_trials=35000 \
      max_minimization_iter=5000 \
      abs_score_convergence_thresh=200.0 \
      backrub_trajectory_stride=35000 \
  -fa_max_dis 9.0 \
  -ex1 -ex2 \
  -beta_nov16 \
  -include_sugars \
  -alternate_3_letter_codes pdb_sugar \
  -maintain_links \
  -ignore_unrecognized_res false \
  -ignore_zero_occupancy false

for pdb in ./*.pdb; do
  [[ -e "$pdb" ]] || continue
  cp -f "$pdb" "$SLURM_SUBMIT_DIR/$SNAPSHOT_DIR/${MUT_TAG}_$(basename "$pdb")"
done

[[ -f score.sc   ]] && cp -f score.sc   "$SLURM_SUBMIT_DIR/$SCORE_DIR/${MUT_TAG}_score.sc"
[[ -f ddG.db3    ]] && cp -f ddG.db3    "$SLURM_SUBMIT_DIR/$DDG_DIR/${MUT_TAG}_ddG.db3"
[[ -f struct.db3 ]] && cp -f struct.db3 "$SLURM_SUBMIT_DIR/$STRUCT_DIR/${MUT_TAG}_struct.db3"

info "Finished Flex ddG for mutation: $MUT_TAG"
