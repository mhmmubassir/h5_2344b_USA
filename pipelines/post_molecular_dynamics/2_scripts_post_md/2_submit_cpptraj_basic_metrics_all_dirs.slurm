#!/usr/bin/env bash
#SBATCH --job-name=postMD_rms_all
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=postMD_rms_all.%j.out
#SBATCH --error=postMD_rms_all.%j.err

bash all_post_md_all_folder.sh
