#!/bin/bash
#SBATCH --job-name="ct_rx_rf15_ha"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=48:00:00
#SBATCH --output=global_relax_%j.out
#SBATCH --error=global_relax_%j.err
#SBATCH --mail-user=mm07305@uga.edu
#SBATCH --mail-type=ALL

module load Rosetta/2022.46.334-intel-2021b

if [[ -z "$1" ]]; then
    echo "Usage: sbatch $0 <protein.pdb>"
    exit 1
fi
protein="$1"
relaxed_protein="relaxed_${protein}"

echo "===================================================="
echo "Running global Cartesian Relax (MPI) on: $protein"
echo "Output file prefix: ${relaxed_protein}"
echo "===================================================="

srun --mpi=pmix_v3 relax.mpi.linuxiccrelease \
  -database /apps/eb/Rosetta/2022.46.334-intel-2021b/database \
  -s "${protein}" \
  -use_input_sc \
  -nstruct 20 \
  -out:file:o "${relaxed_protein}" \
  -relax:cartesian \
  -relax:script cart2.script \
  -relax:min_type lbfgs_armijo_nonmonotone \
  -score:weights ref2015_cart \
  -fa_max_dis 9.0 \
  -constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
  -ideal_sugars \
  -auto_detect_glycan_connections \
  -maintain_links \
  -include_sugars \
  -alternate_3_letter_codes pdb_sugar \
  -ignore_unrecognized_res false \
  -write_pdb_link_records \
  -ex1 -ex2

echo "Done. Relaxed structure(s) saved as ${relaxed_protein}_*.pdb"
