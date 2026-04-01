#!/usr/bin/env bash
#SBATCH --job-name=af3_rmsd
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/af3_rmsd_%j.out
#SBATCH --error=logs/af3_rmsd_%j.err

set -euo pipefail

module load Anaconda3/2024.02-1

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs results

python 3_af3_rmsd_matrix_basic.py \
  --pdb-dir best_pdbs \
  --output-dir results \
  --reference-substring EPI_ISL_18133029
