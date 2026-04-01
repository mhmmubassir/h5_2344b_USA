# MD post-analysis package

This package is organized into two directories:

## 1_mmgbsa
MM/GBSA setup, submission, collection, and CSV combination.

### Run order
1. `1_run_mmgbsa_all_dirs.sh`
2. `4_submit_mmgbsa.slurm` is called by script 1 inside each `mmgbsa/` folder
3. `5_collect_mmgbsa_files.sh`
4. `6_combine_mmgbsa_to_csv.py`
5. Optional quick summary: `7_collect_mmgbsa_delta_totals.sh`

## 2_post_md
General post-MD trajectory analysis.

### Main workflow
1. `1_run_cpptraj_basic_metrics_all_dirs.sh`
2. or submit with `2_submit_cpptraj_basic_metrics_all_dirs.slurm`
3. `3_collect_frame1000_pdbs.py`
4. `4_collect_rmsd_files.py`
5. `5_collect_rmsf_files.py`

### RBS dynamics workflow
6. `6_run_rbs_geometry_all_dirs.sh`
7. `7_compute_rbs_area_all_runs.py`
8. `8_extract_rbs_representative_pdbs.sh`

### Optional structural exploration
9. `9_run_pca_all_dirs_optional.sh`
10. `10_run_glycan_flex_pipeline_optional.sh`
11. `11_collect_glycan_cloud_pdbs_optional.sh`

## What this already covers
- RMSD
- RMSF
- radius of gyration
- SASA
- hydrogen bonds
- frame PDB extraction
- MM/GBSA and decomposition
- RBS geometry distances
- RBS triangle areas
- representative RBS PDB extraction
- optional PCA
- optional glycan-flex clustering/cloud PDB collection
