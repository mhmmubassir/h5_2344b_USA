#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# Extract PDB snapshots from Rosetta struct.db3 files and rename them into a
# more readable format.
#
# Default input directory : outputs/struct_db3
# Default output directory: outputs/extracted_pdbs
#
# Optional environment-variable overrides:
#   ROSETTA_MODULE=Rosetta/2022.46.334-intel-2021b
#   STRUCT_DB_DIR=outputs/struct_db3
#   EXTRACT_OUT_DIR=outputs/extracted_pdbs
#   BACKRUB_STRIDE=35000
# -----------------------------------------------------------------------------

set -euo pipefail
shopt -s nullglob

readonly ROSETTA_MODULE="${ROSETTA_MODULE:-Rosetta/2022.46.334-intel-2021b}"
readonly STRUCT_DB_DIR="${STRUCT_DB_DIR:-outputs/struct_db3}"
readonly EXTRACT_OUT_DIR="${EXTRACT_OUT_DIR:-outputs/extracted_pdbs}"
readonly BACKRUB_STRIDE="${BACKRUB_STRIDE:-35000}"

usage() {
  echo "Usage: STRUCT_DB_DIR=outputs/struct_db3 EXTRACT_OUT_DIR=outputs/extracted_pdbs bash $0"
}

die() {
  echo "❌ $*" >&2
  exit 1
}

info() {
  echo "[$(date '+%F %T')] $*"
}

format_name() {
  local struct_id=$1
  local mutation_tag=$2
  local suffix=$3
  local pose
  local snapshot_number

  case $(( (struct_id - 1) % 3 )) in
    0) pose="backrub" ;;
    1) pose="wt" ;;
    2) pose="mut" ;;
    *) die "Unexpected struct_id mapping for $struct_id" ;;
  esac

  snapshot_number=$(printf "%05d" $(( ((struct_id - 1) / 3 + 1) * BACKRUB_STRIDE )))
  echo "${mutation_tag}_${pose}_${snapshot_number}${suffix}.pdb"
}

[[ -d "$STRUCT_DB_DIR" ]] || die "Input struct-db directory not found: $STRUCT_DB_DIR"
mkdir -p "$EXTRACT_OUT_DIR"

module load "$ROSETTA_MODULE"
SCORE_JD2="$(command -v score_jd2.mpi.linuxiccrelease || true)"
[[ -n "$SCORE_JD2" && -x "$SCORE_JD2" ]] || die "score_jd2.mpi.linuxiccrelease not found after loading $ROSETTA_MODULE"

mapfile -t DB_FILES < <(find "$STRUCT_DB_DIR" -maxdepth 1 -type f -name '*_struct.db3' | sort -V)
(( ${#DB_FILES[@]} > 0 )) || die "No *_struct.db3 files found in $STRUCT_DB_DIR"

info "Found ${#DB_FILES[@]} struct.db3 files"

for db in "${DB_FILES[@]}"; do
  mutation_tag="$(basename "$db" _struct.db3)"
  output_dir="$EXTRACT_OUT_DIR/$mutation_tag"
  mkdir -p "$output_dir"

  info "Extracting: $mutation_tag"

  mpirun -n 1 "$SCORE_JD2" \
    -inout:dbms:database_name "$db" \
    -in:use_database true \
    -include_sugars \
    -alternate_3_letter_codes pdb_sugar \
    -maintain_links \
    -ignore_unrecognized_res false \
    -out:pdb \
    -out:path:pdb "$output_dir" \
    -mute all

  for pdb in "$output_dir"/*_0001*.pdb; do
    [[ -e "$pdb" ]] || continue
    if [[ "$pdb" =~ /([0-9]+)_0001(_low|_last)?\.pdb$ ]]; then
      struct_id="${BASH_REMATCH[1]}"
      suffix="${BASH_REMATCH[2]}"
      new_name="$(format_name "$struct_id" "$mutation_tag" "$suffix")"
      mv -f "$pdb" "$output_dir/$new_name"
    fi
  done
done

info "Done. Extracted PDBs are under: $EXTRACT_OUT_DIR"
