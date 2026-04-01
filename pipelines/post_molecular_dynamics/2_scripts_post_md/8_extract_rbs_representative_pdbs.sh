#!/usr/bin/env bash
# extract_all_rbs_reps.sh
#
# From the current directory (e.g. MD_9dip23vs26), walk all subdirs,
# find *_rbs folders with rbs_area_p1.dat, compute the representative
# frame (closest to mean area), and extract a protomer1+glycan PDB
# from produ.nc with cpptraj.
#
# Output PDBs go to ./rbs_rep_pdbs/ with unique names:
#   <runname>_p1_frame<frame>.pdb
#
# Assumptions:
#   - Each MD run dir has: produ.nc and cplx.parm7 or CPLX.parm7
#   - Protomer 1: residues 1–496
#   - Protomer 1 glycan: residues 1501–1506
#   - rbs_area_p1.dat is in RUN_BASE_rbs/
#
# Adjust the :1-496,1501-1506 mask if needed.

set -euo pipefail

ROOT="$PWD"
OUTDIR="$ROOT/rbs_rep_pdbs"
mkdir -p "$OUTDIR"

log() { printf '%s\n' "$*" >&2; }

command -v cpptraj >/dev/null 2>&1 || { log "Error: cpptraj not found in PATH"; exit 1; }
command -v python >/dev/null 2>&1 || { log "Error: python not found in PATH"; exit 1; }

log "[INFO] Root: $ROOT"
log "[INFO] Output PDBs → $OUTDIR"

# Find all *_rbs directories under ROOT
while IFS= read -r rbsdir; do
  # rbsdir is absolute or relative path to something like:
  #   /.../0.0.9dip_23.amber_renamegly_rbs
  run_dir="$(dirname "$rbsdir")"
  run_base="$(basename "$run_dir")"

  log ""
  log "[INFO] --- Processing run: $run_base ---"
  log "[INFO] run_dir = $run_dir"
  log "[INFO] rbsdir  = $rbsdir"

  # topology
  parm=""
  if [[ -f "$run_dir/cplx.parm7" ]]; then
    parm="$run_dir/cplx.parm7"
  elif [[ -f "$run_dir/CPLX.parm7" ]]; then
    parm="$run_dir/CPLX.parm7"
  else
    log "[WARN]   No cplx.parm7 or CPLX.parm7 in $run_dir — skipping"
    continue
  fi

  # trajectory
  traj="$run_dir/produ.nc"
  if [[ ! -f "$traj" ]]; then
    log "[WARN]   No produ.nc in $run_dir — skipping"
    continue
  fi

  rbs_area_file="$rbsdir/rbs_area_p1.dat"
  if [[ ! -f "$rbs_area_file" ]]; then
    log "[WARN]   No rbs_area_p1.dat in $rbsdir — skipping"
    continue
  fi

  # Use Python (no numpy) to find frame closest to mean area
  frame=$(python <<EOF
import math

path = r"$rbs_area_file"
frames = []
areas = []

with open(path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            fr = int(parts[0])
            ar = float(parts[1])
        except ValueError:
            continue
        frames.append(fr)
        areas.append(ar)

if not frames:
    raise SystemExit("No valid data in " + path)

mean_area = sum(areas) / len(areas)

best_frame = frames[0]
best_diff = abs(areas[0] - mean_area)

for fr, ar in zip(frames, areas):
    diff = abs(ar - mean_area)
    if diff < best_diff:
        best_diff = diff
        best_frame = fr

print(best_frame)
EOF
)

  if [[ -z "$frame" ]]; then
    log "[WARN]   Could not determine representative frame — skipping"
    continue
  fi

  log "[INFO]   Representative frame for $run_base (protomer 1) = $frame"

  # Output PDB name; should be unique by run_base + frame
  out_pdb="$OUTDIR/${run_base}_p1_frame${frame}.pdb"
  log "[INFO]   Will write PDB → $out_pdb"

  # Write cpptraj input in the run_dir
  cpp_in="$run_dir/extract_rbs_rep_p1.in"

  cat > "$cpp_in" <<EOF
parm $parm
trajin $traj

autoimage anchor @CA
rms first :1-1488@N,CA,C

# Keep only protomer 1 (1–496) + its glycan (1501–1506)
strip !(:1-496,1501-1506)

trajout $out_pdb onlyframes $frame
run
quit
EOF

  # Run cpptraj in the run_dir
  (
    cd "$run_dir"
    log "[INFO]   Running cpptraj in $run_dir"
    cpptraj -i "$(basename "$cpp_in")" > extract_rbs_rep_p1.log 2>&1
  )

  if [[ -f "$out_pdb" ]]; then
    log "[INFO]   ✅ Extracted $out_pdb"
  else
    log "[WARN]   cpptraj finished but $out_pdb not found"
  fi

done < <(find "$ROOT" -type d -name '*_rbs')

log ""
log "[INFO] Done. Representative RBS PDBs are in: $OUTDIR"
