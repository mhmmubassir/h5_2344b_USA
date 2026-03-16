#!/bin/bash
#───────────────────────────────────────────────────────────
#           SLURM Job Configuration  
#───────────────────────────────────────────────────────────
#SBATCH --job-name=flex_9dip23_2114
#SBATCH --partition=bahl_p
#SBATCH --constraint=EDR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1     
#SBATCH --mem-per-cpu=5G      
#SBATCH --time=96:00:00
#SBATCH --array=1-2114
#SBATCH --output=logs/flexddg_%A_%a.out
#SBATCH --error=logs/flexddg_%A_%a.err

set -euo pipefail

#───────────────────────────────────────────────────────────
# 0) Directory layout
#───────────────────────────────────────────────────────────
cd "$SLURM_SUBMIT_DIR" || { echo "❌ Cannot cd \$SLURM_SUBMIT_DIR"; exit 1; }
mkdir -p runs structures \
         outputs/{ddg_db3,score_sc,struct_db3} \
         logs

#───────────────────────────────────────────────────────────
# 1) Rosetta
#───────────────────────────────────────────────────────────
module load Rosetta/2022.46.334-intel-2021b

#───────────────────────────────────────────────────────────
# 2) CLI check
#───────────────────────────────────────────────────────────
if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch --array=1-N $0 <pre-minimised_complex.pdb>"
  exit 1
fi
glyco_pdb=$1
[[ -f $glyco_pdb ]] || { echo "❌ PDB '$glyco_pdb' not found"; exit 1; }

#───────────────────────────────────────────────────────────
# 3) Pick resfile for this task
#───────────────────────────────────────────────────────────
mapfile -t res_files < <(ls resfiles/*.resfile | sort -V)
(( ${#res_files[@]} > 0 )) || { echo "❌ No *.resfile in resfiles/"; exit 1; }

tid=${SLURM_ARRAY_TASK_ID:-0}
(( tid>=1 && tid<=${#res_files[@]} )) || { echo "Bad SLURM_ARRAY_TASK_ID=$tid"; exit 1; }

res_file=${res_files[tid-1]}
res_tag=$(basename "$res_file" .resfile)

#───────────────────────────────────────────────────────────
# 4) Scratch working dir
#───────────────────────────────────────────────────────────
run_dir="runs/$res_tag"
mkdir -p "$run_dir"
cp "$glyco_pdb" "$res_file" ddG-backrub.xml "$run_dir/"
cd "$run_dir"

echo -e "\n──────────────────────────────────────────────"
echo "  Mutation tag  : $res_tag"
echo "  SLURM index   : $tid / ${#res_files[@]}"
echo -e "──────────────────────────────────────────────\n"

#───────────────────────────────────────────────────────────
# 5) Flex ddG
#───────────────────────────────────────────────────────────
srun --mpi=pmix_v3 -n "$SLURM_NTASKS" \
     rosetta_scripts.mpi.linuxiccrelease \
     -database /apps/eb/Rosetta/2022.46.334-intel-2021b/database \
     -in:file:s  "$(basename "$glyco_pdb")" \
     -in:file:fullatom \
     -parser:protocol  ddG-backrub.xml \
     -parser:script_vars chainstomove=A \
                         mutate_resfile_relpath="$(basename "$res_file")" \
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

#───────────────────────────────────────────────────────────
# 6) Rename & move outputs
#───────────────────────────────────────────────────────────
# 6a) PDB snapshots
for pdb in *.pdb; do
    cp "$pdb"  "${SLURM_SUBMIT_DIR}/structures/${res_tag}_${pdb}"
done

# 6b) score.sc, ddG.db3, struct.db3
[[ -f score.sc    ]] && mv -f score.sc    "${res_tag}_score.sc"
[[ -f ddG.db3     ]] && mv -f ddG.db3     "${res_tag}_ddG.db3"
[[ -f struct.db3  ]] && mv -f struct.db3  "${res_tag}_struct.db3"

for f in *_score.sc *_ddG.db3 *_struct.db3; do
    case $f in
        *_score.sc)    cp "$f" "${SLURM_SUBMIT_DIR}/outputs/score_sc/$f"    ;;
        *_ddG.db3)     cp "$f" "${SLURM_SUBMIT_DIR}/outputs/ddg_db3/$f"     ;;
        *_struct.db3)  cp "$f" "${SLURM_SUBMIT_DIR}/outputs/struct_db3/$f"  ;;
    esac
done

#───────────────────────────────────────────────────────────
# 7) Copy-then-delete SLURM logs (avoids race w/ SLURM)
#───────────────────────────────────────────────────────────
orig_out="${SLURM_SUBMIT_DIR}/logs/flexddg_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
orig_err="${SLURM_SUBMIT_DIR}/logs/flexddg_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

new_out="${SLURM_SUBMIT_DIR}/logs/${res_tag}.out"
new_err="${SLURM_SUBMIT_DIR}/logs/${res_tag}.err"

cp "$orig_out" "$new_out"  && rm -f "$orig_out"
cp "$orig_err" "$new_err"  && rm -f "$orig_err"

echo "✅  Flex ddG finished for $res_tag"
