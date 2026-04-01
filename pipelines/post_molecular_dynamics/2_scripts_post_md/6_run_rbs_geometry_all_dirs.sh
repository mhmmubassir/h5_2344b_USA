#!/usr/bin/env bash
# rbs_geometry_all_directory.sh — run RBS geometry analysis (distances only)
# Assumes:
#   - Protein: residues 1–1488 (HA trimer, 3 × 496 residues)
#   - Glycans: 1489–1506 (not used here)
#   - Protomer 1:   1–496
#   - Protomer 2: 497–992  (= +496)
#   - Protomer 3: 993–1488 (= +992)
#
# RBS elements (per protomer, protomer 1 numbering):
#   130-loop   ≈ 130–138
#   150-loop   ≈ 147–155  (not used directly here but useful for RMSF later)
#   190-helix  ≈ 187–194
#   220-loop   ≈ 218–225
#   Anchors for triangle: 133 (130-loop), 187 (190-helix), 223 (220-loop)
#
# Usage:
#   bash rbs_geometry_all_directory.sh
#   bash rbs_geometry_all_directory.sh --dry-run

set -euo pipefail

EXCLUDE_DIRS=(y.input_bondfiles y.input_files z.template)
DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

log() { printf '%s\n' "$*" >&2; }

write_and_run_cpptraj() {
  local workdir="$1" parm_path="$2" traj_path="$3"
  local base; base="$(basename "$workdir")"
  local OUTDIR="$workdir/${base}_rbs"
  mkdir -p "$OUTDIR"

  cat > "$OUTDIR/rbs_geometry.in" <<EOF
parm    ${parm_path}
trajin  ${traj_path}
autoimage anchor @CA

###############################################
# RBS geometry metrics (per protomer)
# Protomer 1: residues 1–496
# Protomer 2: residues 497–992  (= protomer 1 + 496)
# Protomer 3: residues 993–1488 (= protomer 1 + 992)
###############################################

# --- Protomer 1 (1–496) ---

# 130-loop vs 220-loop COM distance (mouth width)
distance d_p1_130_220 :130-138 :218-225 out d_p1_130_220.dat

# Triangle sides between anchor CAs (133, 187, 223)
distance d_p1_133_187 :133@CA :187@CA out d_p1_133_187.dat
distance d_p1_187_223 :187@CA :223@CA out d_p1_187_223.dat
distance d_p1_133_223 :133@CA :223@CA out d_p1_133_223.dat

# --- Protomer 2 (497–992 = +496) ---

distance d_p2_130_220 :626-634 :714-721 out d_p2_130_220.dat

distance d_p2_133_187 :629@CA :683@CA out d_p2_133_187.dat
distance d_p2_187_223 :683@CA :719@CA out d_p2_187_223.dat
distance d_p2_133_223 :629@CA :719@CA out d_p2_133_223.dat

# --- Protomer 3 (993–1488 = +992) ---

distance d_p3_130_220 :1122-1130 :1210-1217 out d_p3_130_220.dat

distance d_p3_133_187 :1125@CA :1179@CA out d_p3_133_187.dat
distance d_p3_187_223 :1179@CA :1215@CA out d_p3_187_223.dat
distance d_p3_133_223 :1125@CA :1215@CA out d_p3_133_223.dat

# --- Optional: RBS loop RMSF (uncomment if/when needed) ---
# atomicfluct byres out rmsf_rbs_130loop.dat  :130-138,626-634,1122-1130@N,CA,C
# atomicfluct byres out rmsf_rbs_150loop.dat  :147-155,643-651,1139-1147@N,CA,C
# atomicfluct byres out rmsf_rbs_190helix.dat :187-194,683-690,1179-1186@N,CA,C
# atomicfluct byres out rmsf_rbs_220loop.dat  :218-225,714-721,1210-1217@N,CA,C

run
quit
EOF

  if $DRY_RUN; then
    log "[DRY] Would run cpptraj in: $OUTDIR"
    return 0
  fi

  (
    cd "$OUTDIR"
    log "[INFO] ($base) Running cpptraj RBS geometry…"
    cpptraj -i rbs_geometry.in > rbs_geometry.log 2>&1
  )

  log "✅ ($base) Done → $OUTDIR"
  tail -n 10 "$OUTDIR/rbs_geometry.log" || true
}

# ---- main -----------------------------------------------------------

command -v cpptraj >/dev/null 2>&1 || { log "Error: cpptraj not found in PATH"; exit 1; }

shopt -s nullglob
found_any=false

for d in "$PWD"/*/; do
  found_any=true
  dir="${d%/}"
  base="$(basename "$dir")"

  # skip excluded
  skip=false
  for ex in "${EXCLUDE_DIRS[@]}"; do
    [[ "$base" == "$ex" ]] && { skip=true; break; }
  done
  $skip && { log "[SKIP] $base — in exclude list"; continue; }

  # find parm and traj
  parm=""
  [[ -f "$dir/cplx.parm7" ]] && parm="$dir/cplx.parm7"
  [[ -z "$parm" && -f "$dir/CPLX.parm7" ]] && parm="$dir/CPLX.parm7"

  traj=""
  [[ -f "$dir/produ.nc" ]] && traj="$dir/produ.nc"

  if [[ -z "$parm" || -z "$traj" ]]; then
    $DRY_RUN && log "[DRY][SKIP] $base — missing cplx.parm7/CPLX.parm7 or produ.nc"
    [[ -z "$parm" ]] && log "[SKIP] $dir — missing cplx.parm7 (or CPLX.parm7)"
    [[ -z "$traj" ]] && log "[SKIP] $dir — missing produ.nc"
    continue
  fi

  log "[INFO] ($base) Using parm=$parm traj=$traj"
  write_and_run_cpptraj "$dir" "$parm" "$traj"
done

if ! $found_any; then
  log "No immediate subdirectories found. Nothing to do."
fi
