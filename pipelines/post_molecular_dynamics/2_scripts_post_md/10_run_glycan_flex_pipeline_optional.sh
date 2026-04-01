#!/usr/bin/env bash
# 30_all_glycan_flex.sh (PDB-only, robust version)
# Run from MD_9dip23vs26/
#  A: per-run RBS + nearest-glycan prep → multi-MODEL PDB
#  B: per-variant×linkage clustering + PCA (from those PDBs)
#  C: representative cluster-1 PDBs for PyMOL

set -u
set -o pipefail

# ----------------------------------------------------------------------
# Check cpptraj
# ----------------------------------------------------------------------
if ! command -v cpptraj >/dev/null 2>&1; then
  echo "[ERROR] cpptraj not found in PATH. Did you 'ml amber' first?" >&2
  exit 1
fi

ROOT="$PWD"
OUT="$ROOT/glycan_flex"
PREP="$OUT/prep_per_run"
VAROUT="$OUT/by_variant"
REP="$OUT/representatives"

mkdir -p "$PREP" "$VAROUT" "$REP"

MAPFILE="$OUT/variant_run_list.tsv"
: > "$MAPFILE"

echo "[INFO] Root:   $ROOT"
echo "[INFO] Output: $OUT"

# ----------------------------------------------------------------------
# Masks & parameters
# ----------------------------------------------------------------------
RBS_VIEW_RANGE="90-230"

# Glycan residue ranges for F/G/H (protomer 1)
GLY_F_RANGE="1489-1494"   # core 1489-1493 + SA 1494
GLY_G_RANGE="1495-1500"   # core 1495-1499 + SA 1500
GLY_H_RANGE="1501-1506"   # core 1501-1505 + SA 1506

CLUSTER_RMS_CUTOFF="1.5"  # Å; adjust if needed

# For clustering/PCA: use all heavy atoms in stripped RBS+glycan
RMS_MASK="!@H="

# ======================================================================
# STAGE A – per-run prep: detect closest glycan + strip to RBS+glycan
# ======================================================================

echo "[INFO] Stage A: per-run prep with nearest-glycan detection"

mapfile -t RBS_DIRS < <(find "$ROOT" -maxdepth 6 -type d -name '*_rbs' | sort)

if [[ ${#RBS_DIRS[@]} -eq 0 ]]; then
  echo "[ERROR] No *_rbs directories found. Are you in MD_9dip23vs26?" >&2
  exit 1
fi

echo "[INFO] Found ${#RBS_DIRS[@]} *_rbs directories."

runs_processed=0

for rdir in "${RBS_DIRS[@]}"; do
  run_dir="$(dirname "$rdir")"     # e.g. .../0.0.9dip_23.amber_renamegly
  run_name="$(basename "$run_dir")"

  # Strip first two dot-separated fields: "0.0." / "1.2." etc.
  core="${run_name#*.*.}"
  # core examples:
  #   0.0.9dip_23.amber_renamegly      -> 9dip_23.amber_renamegly
  #   0.0.2021_ref_23                  -> 2021_ref_23
  #   1.0.final_MUT_..._V3_26_amber    -> final_MUT_..._V3_26_amber

  if [[ "$core" == *_23* ]]; then
    linkage="23"
    variant="${core%%_23*}"
  elif [[ "$core" == *_26* ]]; then
    linkage="26"
    variant="${core%%_26*}"
  else
    linkage="unk"
    variant="$core"
  fi

  tag="${variant}_${linkage}"

  # Parm/trajectory
  parm=""
  if [[ -f "$run_dir/CPLX.parm7" ]]; then
    parm="$run_dir/CPLX.parm7"
  elif [[ -f "$run_dir/cplx.parm7" ]]; then
    parm="$run_dir/cplx.parm7"
  else
    echo "[WARN] $run_dir: no CPLX.parm7/cplx.parm7; skipping."
    continue
  fi

  traj="$run_dir/produ.nc"
  if [[ ! -f "$traj" ]]; then
    echo "[WARN] $run_dir: no produ.nc; skipping."
    continue
  fi

  # ---- A1. Detect which sialic acid (F/G/H) is closest to RBS (frame 1) ----
  dF_dat="$PREP/${run_name}_dF.dat"
  dG_dat="$PREP/${run_name}_dG.dat"
  dH_dat="$PREP/${run_name}_dH.dat"
  detect_in="$PREP/${run_name}_detect_gly.in"
  detect_log="$PREP/${run_name}_detect_gly.log"

  cat > "$detect_in" <<EOF
parm $parm
trajin $traj 1 1 1

autoimage anchor @CA
center :${RBS_VIEW_RANGE}@CA

distance dF :${RBS_VIEW_RANGE}@CA :1494&!@H= out $dF_dat
distance dG :${RBS_VIEW_RANGE}@CA :1500&!@H= out $dG_dat
distance dH :${RBS_VIEW_RANGE}@CA :1506&!@H= out $dH_dat

run
quit
EOF

  if ! cpptraj -i "$detect_in" > "$detect_log" 2>&1; then
    echo "[WARN] $run_name: cpptraj detect_gly failed, see $detect_log; skipping run."
    continue
  fi

  dF=$(awk 'NF>1 && $1 !~ /^#/ {print $2; exit}' "$dF_dat" 2>/dev/null || echo "9999")
  dG=$(awk 'NF>1 && $1 !~ /^#/ {print $2; exit}' "$dG_dat" 2>/dev/null || echo "9999")
  dH=$(awk 'NF>1 && $1 !~ /^#/ {print $2; exit}' "$dH_dat" 2>/dev/null || echo "9999")

  closest=$(awk -v a="$dF" -v b="$dG" -v c="$dH" 'BEGIN{
      min=a; g="F";
      if (b < min) {min=b; g="G"}
      if (c < min) {min=c; g="H"}
      print g
  }')

  case "$closest" in
    F) GLY_RANGE="$GLY_F_RANGE" ;;
    G) GLY_RANGE="$GLY_G_RANGE" ;;
    H) GLY_RANGE="$GLY_H_RANGE" ;;
    *) echo "[WARN] $run_name: could not determine glycan (F/G/H); skipping run."; continue ;;
  esac

  echo "[INFO] $run_name → tag=$tag, closest glycan=$closest (resi $GLY_RANGE)"

  # ---- A2. Strip to RBS+that glycan, align, write multi-MODEL PDB ---------
  pdb_out="$PREP/${run_name}_rbs_gly.pdb"
  prep_in="$PREP/${run_name}_prep.in"
  prep_log="$PREP/${run_name}_prep.log"

  cat > "$prep_in" <<EOF
parm $parm
trajin $traj

autoimage anchor @CA

# keep RBS + selected glycan only
strip !(:${RBS_VIEW_RANGE},:${GLY_RANGE})

center :${RBS_VIEW_RANGE}@CA
rms first :${RBS_VIEW_RANGE}@N,CA,C

trajout $pdb_out pdb
run
quit
EOF

  if ! cpptraj -i "$prep_in" > "$prep_log" 2>&1; then
    echo "[WARN] $run_name: cpptraj prep failed, see $prep_log; skipping run."
    continue
  fi

  if [[ ! -s "$pdb_out" ]]; then
    echo "[WARN] $run_name: $pdb_out not created; skipping run."
    continue
  fi

  echo -e "${tag}\t${pdb_out}" >> "$MAPFILE"
  ((runs_processed++))
done

echo "[INFO] Stage A complete. Runs processed: $runs_processed"
echo "[INFO] Mapfile: $MAPFILE"

if [[ $runs_processed -eq 0 ]]; then
  echo "[ERROR] No runs were successfully processed in Stage A." >&2
  exit 1
fi

# ======================================================================
# STAGE B – per-variant clustering + PCA (from PDBs)
# ======================================================================

echo "[INFO] Stage B: clustering + PCA per variant×linkage"

mapfile -t TAGS < <(awk '{print $1}' "$MAPFILE" | sort -u)

for tag in "${TAGS[@]}"; do
  mapfile -t LINES < <(awk -v t="$tag" '$1==t{print $0}' "$MAPFILE")
  if [[ ${#LINES[@]} -eq 0 ]]; then
    continue
  fi

  var_dir="$VAROUT/$tag"
  mkdir -p "$var_dir"

  cpp_in="$var_dir/${tag}_cluster_pca.in"
  cpp_log="$var_dir/${tag}_cluster_pca.log"
  ensemble_pdb="$var_dir/${tag}_ensemble.pdb"
  rep_pdb="$var_dir/${tag}_cluster_rep.pdb"
  cluster_dat="$var_dir/${tag}_cluster.dat"
  cluster_sum="$var_dir/${tag}_cluster_summary.dat"
  pca_dat="$var_dir/${tag}_pca.dat"

  # use the first PDB for this tag as parm
  first_pdb=$(echo "${LINES[0]}" | awk '{print $2}')

  {
    echo "parm $first_pdb"
    for line in "${LINES[@]}"; do
      pdb_file=$(echo "$line" | awk '{print $2}')
      echo "trajin $pdb_file"
    done

    echo "autoimage anchor @CA"
    echo "center :${RBS_VIEW_RANGE}@CA"
    echo "rms first ${RMS_MASK}"

    echo "createcrd rbs_gly_crd"

    echo "trajout $ensemble_pdb pdb"

    echo "cluster crd rbs_gly_crd out $cluster_dat summary $cluster_sum \\"
    echo "   repout $rep_pdb repfmt pdb \\"
    echo "   averagelinkage epsilon ${CLUSTER_RMS_CUTOFF} \\"
    echo "   rms ${RMS_MASK}"

    echo "principal name ${tag}_pca crd rbs_gly_crd out $pca_dat"

    echo "run"
    echo "quit"
  } > "$cpp_in"

  echo "[INFO]   $tag: cpptraj clustering/PCA"
  if ! cpptraj -i "$cpp_in" > "$cpp_log" 2>&1; then
    echo "[WARN]   $tag: cpptraj clustering/PCA failed, see $cpp_log; skipping this tag."
    continue
  fi
done

echo "[INFO] Stage B complete."

# ======================================================================
# STAGE C – extract hero representatives (cluster 1) for PyMOL
# ======================================================================

echo "[INFO] Stage C: extracting representative PDBs (cluster 1) for each tag"

for tag in "${TAGS[@]}"; do
  var_dir="$VAROUT/$tag"
  rep_pdb="$var_dir/${tag}_cluster_rep.pdb"

  if [[ ! -f "$rep_pdb" ]]; then
    continue
  fi

  cp "$rep_pdb" "$REP/${tag}_cluster_reps.pdb"

  out1="$REP/${tag}_cluster1.pdb"
  awk '
    /^MODEL/{m++}
    m==1{print}
    /^ENDMDL/ && m==1{exit}
  ' "$rep_pdb" > "$out1"

  echo "[INFO]   → $out1"
done

echo "[INFO] All done."
echo "[INFO] Key outputs:"
echo "  $PREP/                       : per-run stripped RBS+nearest-gly PDBs"
echo "  $VAROUT/<tag>/ensemble       : combined ensembles (PDB)"
echo "  $VAROUT/<tag>/_cluster_*     : cluster populations + medoids + PCA"
echo "  $REP/<tag>_cluster1.pdb      : hero structures for PyMOL"
