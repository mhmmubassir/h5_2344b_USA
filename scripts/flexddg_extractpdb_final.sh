#!/usr/bin/env bash
# extract_all_db3.sh  –  organised pose dump using the MPI score_jd2
# Run inside newtry_9dip23_again

set -euo pipefail

### CONFIG (edit only if these ever change) ############################
ROS_VER="2022.46.334-intel-2021b"   # Rosetta module you used
BACKRUB_STRIDE=35000                # backrub_trajectory_stride in Flex ddG
DB_DIR="ddg_outputs"                # where *_struct.db3 live
OUT_ROOT="extracted_pdbs"           # output root
########################################################################

echo -e "\n>>> Loading Rosetta module"
module load Rosetta/$ROS_VER

echo ">>> Locating MPI score_jd2 binary"
SCORE_JD2=$(command -v score_jd2.mpi.linuxiccrelease) || true
[[ -x $SCORE_JD2 ]] || { echo "   ❌ MPI score_jd2 not found"; exit 1; }
echo "   ✔ $SCORE_JD2"

echo ">>> Scanning $DB_DIR for databases"
mapfile -t DBLIST < <(ls "$DB_DIR"/*_struct.db3 2>/dev/null || true)
[[ ${#DBLIST[@]} -gt 0 ]] || { echo "   ❌ No *_struct.db3 files found"; exit 1; }
echo "   ✔ Found ${#DBLIST[@]} databases"

mkdir -p "$OUT_ROOT"

# helper → descriptive filename
fname() {                 # $1 struct_id, $2 muttag, $3 suffix
  case $(( ( $1 - 1 ) % 3 )) in
     0) pose=backrub ;;
     1) pose=wt ;;
     2) pose=mut ;;
  esac
  num=$(printf "%05d" $(( ( ($1-1)/3 +1 ) * BACKRUB_STRIDE )))
  echo "${2}_${pose}_${num}${3}.pdb"
}

echo ">>> Extracting PDBs"
for db in "${DBLIST[@]}"; do
  muttag=$(basename "$db" _struct.db3)
  outdir="$OUT_ROOT/$muttag"; mkdir -p "$outdir"
  echo "   → $muttag"

  I_MPI_PMI_LIBRARY=mpich mpirun -n 1 "$SCORE_JD2" \
      -inout:dbms:database_name "$db" \
      -in:use_database true \
      -include_sugars \
      -alternate_3_letter_codes pdb_sugar \
      -maintain_links \
      -ignore_unrecognized_res false \
      -out:pdb -out:path:pdb "$outdir" \
      -mute all

  for f in "$outdir"/*_0001*.pdb; do
      [[ $f =~ /([0-9]+)_0001(_low|_last)?\.pdb$ ]] || continue
      id=${BASH_REMATCH[1]}  suf=${BASH_REMATCH[2]}
      mv -f "$f" "$outdir/$(fname "$id" "$muttag" "$suf")"
  done
done

echo -e "\n✅  Done – PDBs organised under  $OUT_ROOT/<mutation>/\n"
