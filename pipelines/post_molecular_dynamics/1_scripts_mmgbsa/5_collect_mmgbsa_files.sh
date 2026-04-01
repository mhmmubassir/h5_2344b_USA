#!/usr/bin/env bash
set -euo pipefail

# Run this script FROM INSIDE MD_9dip23vs26
ROOT="."                         # current dir = MD_9dip23vs26
DEST="$ROOT/mmgbsa_collect"
mkdir -p "$DEST"

echo "Collecting MMGBSA files under $ROOT ..."
echo "Destination: $DEST"
echo

# Find all candidate mmgbsa dirs (recursively, all sub/sub-sub dirs)
mapfile -t MMGBSA_DIRS < <(
  find "$ROOT" -type d \( -name "mmgbsa" -o -name "new_mmgbsa" -o -name "new_mmgbsalab" \) | sort
)

if [ ${#MMGBSA_DIRS[@]} -eq 0 ]; then
  echo "No mmgbsa/new_mmgbsa/new_mmgbsalab directories found under $ROOT"
  exit 1
fi

for src in "${MMGBSA_DIRS[@]}"; do
  base=$(basename "$src")

  # If this is a parent mmgbsa and it contains a new_* child, skip the parent
  if [[ "$base" == "mmgbsa" ]]; then
    if [ -d "$src/new_mmgbsa" ] || [ -d "$src/new_mmgbsalab" ]; then
      echo "Skipping parent mmgbsa (new_* exists): $src"
      continue
    fi
  fi

  # Skip OLD 9dip26 2021ref runs (0.0/0.1/0.2.2021_ref_26)
  if [[ "$src" == *"/9dip26/mdrun_500ns_2021ref/0."* ]]; then
    echo "Skipping old 9dip26 2021ref dir: $src"
    continue
  fi

  # Make a label based on the relative path
  rel="${src#./}"               # drop leading "./"
  label="${rel//\//_}"          # replace / with _

  outdir="$DEST/$label"
  mkdir -p "$outdir"

  echo ">>> $src"
  echo "    -> $outdir"

  # Copy only the needed .dat files (no CSV)
  for f in mmgbsa_gbneck2.dat decom_gbneck2.dat; do
    if [ -f "$src/$f" ]; then
      cp -v "$src/$f" "$outdir/${label}_$f"
    fi
  done

  echo
done

echo "Done."
