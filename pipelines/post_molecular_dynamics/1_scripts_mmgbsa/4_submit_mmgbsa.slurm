#!/bin/bash
#SBATCH --job-name=MMGBSA
#SBATCH --partition=batch
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
#SBATCH --output=MMGBSA.%j.out
#SBATCH --error=MMGBSA.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mm07305@uga.edu

module load Amber/24.3-foss-2022a-AmberTools-24.10-CUDA-12.1.1
source ${AMBERHOME}/amber.sh
cd   $SLURM_SUBMIT_DIR

mpirun -np $SLURM_NTASKS  MMPBSA.py.MPI -O \
      -i  mmgbsa.in \
      -cp CPLX.parm7 \
      -rp receptor.parm7 \
      -lp ligand.parm7 \
      -y  produ.nc \
      -o  mmgbsa_gbneck2.dat \
      -do decom_gbneck2.dat \
      -deo per_frame_breakdown.csv
