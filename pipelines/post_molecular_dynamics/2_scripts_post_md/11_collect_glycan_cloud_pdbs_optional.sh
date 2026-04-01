#!/usr/bin/env bash
# 31_collect_cloud_pdbs.sh
# Run from: MD_9dip23vs26
# Collect all per-run RBS+nearest-glycan PDBs into one folder for PyMOL “cloud” views.

set -euo pipefail

ROOT="$(pwd)"
SRC="${ROOT}/glycan_flex/prep_per_run"
DEST="${ROOT}/glycan_flex/cloud_pdbs"

echo "[INFO] Root:   ${ROOT}"
echo "[INFO] Source: ${SRC}"
echo "[INFO] Dest:   ${DEST}"

mkdir -p "${DEST}"

shopt -s nullglob
count=0
for pdb in "${SRC}"/*_rbs_gly.pdb; do
    base="$(basename "${pdb}")"
    echo "[COPY] ${base}"
    cp -f "${pdb}" "${DEST}/${base}"
    ((count++))
done

if (( count == 0 )); then
    echo "[WARN] No *_rbs_gly.pdb found in ${SRC}"
else
    echo "[INFO] Copied ${count} PDBs into ${DEST}"
fi
