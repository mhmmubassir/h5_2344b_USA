#!/usr/bin/env bash
# all_postMDrms_analysis_all_directory.sh — strict one-level batch for cpptraj
# Only runs in immediate subdirs that contain:
#   cplx.parm7 (or CPLX.parm7)  AND  produ.nc
# Skips: y.input_bondfiles, y.input_files, z.template
# Usage:
#   bash all_postMDrms_analysis_all_directory.sh         # (on compute node via sbatch)
#   bash all_postMDrms_analysis_all_directory.sh --dry-run
set -euo pipefail

EXCLUDE_DIRS=(y.input_bondfiles y.input_files z.template)
DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

log() { printf '%s\n' "$*" >&2; }

write_and_run_cpptraj() {
  local workdir="$1" parm_path="$2" traj_path="$3"
  local base; base="$(basename "$workdir")"
  local OUTDIR="$workdir/${base}_analysis"
  mkdir -p "$OUTDIR"

  cat > "$OUTDIR/analysis_noPCA.in" <<EOF
###############################################
#  analysis_noPCA.in (auto-generated)
#  SYSTEM: protein (1–1488) + three 6-sugar glycans (1489–1506)
###############################################

parm    ${parm_path}
trajin  ${traj_path}
autoimage anchor @CA

# Metrics
rms first out rmsd_protein.dat            :1-1488@N,CA,C
rms first out rmsd_glycore_F.dat          :1489-1493&!@H=
rms first out rmsd_glycore_G.dat          :1495-1499&!@H=
rms first out rmsd_glycore_H.dat          :1501-1505&!@H=
rms first out rmsd_sa_F.dat               :1494&!@H=
rms first out rmsd_sa_G.dat               :1500&!@H=
rms first out rmsd_sa_H.dat               :1506&!@H=
rms first out rmsd_all_sa.dat             :1494,1500,1506&!@H=
rms first out rmsd_all_glycan.dat         :1489-1506&!@H=
rms first out rmsd_complex.dat            (:1-1488@N,CA,C)|(:1489-1506&!@H=)

atomicfluct out rmsf_protein.dat          byres :1-1488@N,CA,C
atomicfluct out rmsf_glycore_F.dat        byres :1489-1493&!@H=
atomicfluct out rmsf_glycore_G.dat        byres :1495-1499&!@H=
atomicfluct out rmsf_glycore_H.dat        byres :1501-1505&!@H=
atomicfluct out rmsf_sa_F.dat             byres :1494&!@H=
atomicfluct out rmsf_sa_G.dat             byres :1500&!@H=
atomicfluct out rmsf_sa_H.dat             byres :1506&!@H=
atomicfluct out rmsf_complex.dat          byres (:1-1488@N,CA,C)|(:1489-1506&!@H=)

radgyr out radgyr_protein.dat             :1-1488
radgyr out radgyr_glycore_F.dat           :1489-1493
radgyr out radgyr_glycore_G.dat           :1495-1499
radgyr out radgyr_glycore_H.dat           :1501-1505
radgyr out radgyr_sa_F.dat                :1494
radgyr out radgyr_sa_G.dat                :1500
radgyr out radgyr_sa_H.dat                :1506
radgyr out radgyr_complex.dat             :1-1506

surf :1-1488                  out sasa_protein.dat
surf :1489-1505               out sasa_glycore_all.dat
surf :1494,1500,1506          out sasa_all_sa.dat
surf :1-1506                  out sasa_complex.dat

hbond out hbonds_protein_Fgly.dat   series donormask :1-1488    acceptormask :1489-1494 avgout hbonds_protein_Fgly_avg.dat
hbond out hbonds_protein_Ggly.dat   series donormask :1-1488    acceptormask :1495-1500 avgout hbonds_protein_Ggly_avg.dat
hbond out hbonds_protein_Hgly.dat   series donormask :1-1488    acceptormask :1501-1506 avgout hbonds_protein_Hgly_avg.dat
hbond out hbonds_protein_allgly.dat series donormask :1-1488    acceptormask :1489-1506 avgout hbonds_protein_allgly_avg.dat
hbond out hbonds_Fcore_Fsa.dat      series donormask :1489-1493 acceptormask :1494      avgout hbonds_Fcore_Fsa_avg.dat
hbond out hbonds_Gcore_Gsa.dat      series donormask :1495-1499 acceptormask :1500      avgout hbonds_Gcore_Gsa_avg.dat
hbond out hbonds_Hcore_Hsa.dat      series donormask :1501-1505 acceptormask :1506      avgout hbonds_Hcore_Hsa_avg.dat

# Frame outputs (keep light; all_frames.pdb is commented out)
# trajout all_frames.pdb
trajout frame001.pdb  onlyframes 1
trajout frame250.pdb  onlyframes 250
trajout frame500.pdb  onlyframes 500
trajout frame750.pdb  onlyframes 750
trajout frame1000.pdb onlyframes 1000

run
quit
EOF

  if $DRY_RUN; then
    log "[DRY] Would run cpptraj in: $OUTDIR"
    return 0
  fi

  (
    cd "$OUTDIR"
    log "[INFO] ($(basename "$workdir")) Running cpptraj…"
    cpptraj -i analysis_noPCA.in > analysis.log 2>&1
  )

  log "✅ ($(basename "$workdir")) Done → $OUTDIR"
  tail -n 12 "$OUTDIR/analysis.log" || true
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

  # require exact filenames directly inside the subdir
  parm=""
  [[ -f "$dir/cplx.parm7" ]] && parm="$dir/cplx.parm7"
  [[ -z "$parm" && -f "$dir/CPLX.parm7" ]] && parm="$dir/CPLX.parm7"

  traj=""
  [[ -f "$dir/produ.nc" ]] && traj="$dir/produ.nc"

  if [[ -z "$parm" || -z "$traj" ]]; then
    $DRY_RUN && log "[DRY][SKIP] $base — missing required files (need cplx.parm7 or CPLX.parm7 AND produ.nc)"
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
