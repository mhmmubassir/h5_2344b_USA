# AF3 basic protocol

This is the simplified AF3 workflow for the version you actually used:

1. Run AF3 prediction jobs in blocks
2. Compute pairwise whole-HA and RBS RMSD matrices
3. Compute RBS opening / loop geometry features
4. Plot AF3 RBS geometry figures

## Script order

1. `1_af3_submit_blocks_basic.sh`
2. `2_af3_run_basic.slurm`
3. `3_af3_rmsd_matrix_basic.py`
4. `4_af3_rmsd_matrix_basic.slurm`
5. `5_af3_rbs_open_basic.py`
6. `6_af3_plot_rbs_geometry_basic.py`


### AF3 prediction
- `input.lst`
- `jsons/`

### RMSD
- `best_pdbs/`

### RBS opening
- `best_pdbs/`
- `compare_numbering.csv`

### plotting
- `rbs_opening_features.csv`
- `top_variants_membership.csv`
- `tips_lineageA_lineageB_mapping.csv`

