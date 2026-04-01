#!/usr/bin/env bash
# pca_all_dirs.sh — batch PCA (cpptraj + Python density plot + frame extract)
# - One-level scan of subdirs; requires cplx.parm7 (or CPLX.parm7) AND produ.nc
# - Outputs go to <subdir>/<subdir>_analysis/ (same folder as RMSD script)
# Usage:
#   bash pca_all_dirs.sh            # run
#   bash pca_all_dirs.sh --dry-run  # list what would run
set -euo pipefail

EXCLUDE_DIRS=(y.input_bondfiles y.input_files z.template)
DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

# Optional env overrides (kept for parity with your single-run script)
MASK="${MASK:-@N,C,CA}"
CRDFRAMES_RANGE="${CRDFRAMES_RANGE:-1,1000}"
OUT_PREFIX="${OUT_PREFIX:-fig4.2}"

log() { printf '%s\n' "$*" >&2; }

run_pca_one_dir() {
  local workdir="$1" parm="$2" traj="$3"
  local base; base="$(basename "$workdir")"
  local OUTDIR="$workdir/${base}_analysis"
  mkdir -p "$OUTDIR"

  if $DRY_RUN; then
    log "[DRY] Would run PCA in: $OUTDIR"
    log "[DRY] Using parm=$parm traj=$traj"
    return 0
  fi

  log "[INFO] ($base) [1/3] cpptraj PCA/projection → evecs.dat, pca.dat, normalmodes.nmd, modes.prmtop, modes01.nc"
  cpptraj <<EOF
parm ${parm}
trajin ${traj}
autoimage anchor @CA

rms first ${MASK}
average ${MASK} crdset trajaverage
createcrd trajectory
run

crdaction trajectory rms ref trajaverage ${MASK}
crdaction trajectory matrix covar name covmatrix ${MASK}

runanalysis diagmatrix covmatrix out ${OUTDIR}/evecs.dat vecs 10 name eigenvectors \
  nmwiz nmwizvecs 10 nmwizfile ${OUTDIR}/normalmodes.nmd nmwizmask ${MASK}

crdaction trajectory projection modes eigenvectors beg 1 end 10 ${MASK} out \
  ${OUTDIR}/pca.dat crdframes ${CRDFRAMES_RANGE}

clear all
readdata ${OUTDIR}/evecs.dat name eigenvectors
parm ${parm}
parmwrite out ${OUTDIR}/modes.prmtop
runanalysis modes name eigenvectors trajout ${OUTDIR}/modes01.nc pcmin -150 pcmax 150 \
  tmode 1 trajoutmask ${MASK} trajoutfmt netcdf
quit
EOF

  log "[INFO] ($base) [2/3] Python density plot + most-dense frame → PNG + most_dense_frame.txt"
  (
    cd "$OUTDIR"
    python - <<'PYEOF'
import numpy as np
import matplotlib.pyplot as plt

# Try SciPy KDE; if missing, fall back to 2D hist density
try:
    from scipy.stats import gaussian_kde
    HAVE_KDE = True
except Exception:
    HAVE_KDE = False

# Try seaborn styling if available (optional)
try:
    import seaborn as sns
    sns.set(style='whitegrid', font_scale=1.4)
except Exception:
    pass

data = np.loadtxt('pca.dat')
frame_idx = data[:,0].astype(int)
pc1 = data[:,1]
pc2 = data[:,2]

if HAVE_KDE:
    xy = np.vstack([pc1, pc2])
    z = gaussian_kde(xy)(xy)
else:
    # 2D histogram fallback → map bin counts to each point
    bins = 60
    H, xedges, yedges = np.histogram2d(pc1, pc2, bins=bins)
    # index each point's bin
    xbin = np.searchsorted(xedges, pc1, side='right') - 1
    ybin = np.searchsorted(yedges, pc2, side='right') - 1
    xbin = np.clip(xbin, 0, H.shape[0]-1)
    ybin = np.clip(ybin, 0, H.shape[1]-1)
    z = H[xbin, ybin].astype(float)

imax = int(np.argmax(z))
dense_frame = int(frame_idx[imax])
dense_pc1 = float(pc1[imax])
dense_pc2 = float(pc2[imax])

plt.figure(figsize=(8,6))
sc = plt.scatter(pc1, pc2, c=z, s=50, edgecolor='k', cmap='viridis', alpha=0.75)
plt.scatter(dense_pc1, dense_pc2, color='red', edgecolor='k', s=120, marker='*',
            label=f'Most Dense (Frame {dense_frame})')
cbar = plt.colorbar(sc)
cbar.set_label('Density' + (' (KDE)' if HAVE_KDE else ' (hist)'))

plt.xlabel('Principal Component 1 (PC1)', fontsize=18)
plt.ylabel('Principal Component 2 (PC2)', fontsize=18)
plt.title('PCA Projection with Density Coloring', fontsize=20, weight='bold')
plt.legend(fontsize=12)
try:
    import seaborn as sns
    sns.despine()
except Exception:
    pass
plt.tight_layout()
plt.savefig('fig4.2_pca_projection_density_colored_highlighted_with_legend.png', dpi=300)

with open('most_dense_frame.txt','w') as f:
    f.write(str(dense_frame)+"\n")

print(f"Most Dense Point: PC1={dense_pc1:.4f}, PC2={dense_pc2:.4f}, Frame={dense_frame}")
PYEOF
  )

  local DENSE_FRAME
  DENSE_FRAME="$(tr -d '[:space:]' < "$OUTDIR/most_dense_frame.txt")"
  local OUT_PDB="${OUTDIR}/${OUT_PREFIX}.most_dense_frame${DENSE_FRAME}.pdb"

  log "[INFO] ($base) [3/3] cpptraj: extract frame ${DENSE_FRAME} → $(basename "$OUT_PDB")"
  cpptraj <<EOF
parm ${parm}
trajin ${traj}
autoimage
trajout ${OUT_PDB} onlyframes ${DENSE_FRAME}
run
quit
EOF

  log "✅ ($base) Done → $OUTDIR"
  log "   • Plot:    ${OUTDIR}/fig4.2_pca_projection_density_colored_highlighted_with_legend.png"
  log "   • Dense:   frame ${DENSE_FRAME} (saved in ${OUTDIR}/most_dense_frame.txt)"
  log "   • PDB:     ${OUT_PDB}"
}

# ---- main: iterate one level of subdirs --------------------------------------

command -v cpptraj >/dev/null 2>&1 || { log "Error: cpptraj not found in PATH"; exit 1; }

shopt -s nullglob
found_any=false
for d in "$PWD"/*/; do
  found_any=true
  dir="${d%/}"; base="$(basename "$dir")"

  # skip excluded
  skip=false
  for ex in "${EXCLUDE_DIRS[@]}"; do
    [[ "$base" == "$ex" ]] && { skip=true; break; }
  done
  $skip && { log "[SKIP] $base — in exclude list"; continue; }

  # require exact files in this subdir
  parm=""
  [[ -f "$dir/cplx.parm7" ]] && parm="$dir/cplx.parm7"
  [[ -z "$parm" && -f "$dir/CPLX.parm7" ]] && parm="$dir/CPLX.parm7"
  traj=""
  [[ -f "$dir/produ.nc" ]] && traj="$dir/produ.nc"

  if [[ -z "$parm" || -z "$traj" ]]; then
    $DRY_RUN && log "[DRY][SKIP] $base — missing required files"
    [[ -z "$parm" ]] && log "[SKIP] $dir — missing cplx.parm7 (or CPLX.parm7)"
    [[ -z "$traj" ]] && log "[SKIP] $dir — missing produ.nc"
    continue
  fi

  log "[INFO] ($base) Using parm=$parm traj=$traj"
  run_pca_one_dir "$dir" "$parm" "$traj"
done

if ! $found_any; then
  log "No immediate subdirectories found. Nothing to do."
fi
