# Flex ddG protocol layout

## Suggested repository structure

```text
repo/
├── protocols/
│   └── ddG-backrub.xml
├── inputs/
│   ├── resfiles/
│   └── pdbs/
├── scripts/
│   ├── run_flexddg_array.slurm
│   ├── extract_structures.sh
│   ├── summarize_flexddg.py
│   └── rosetta_to_amber.py
├── outputs/
│   ├── ddg_db3/
│   ├── score_sc/
│   ├── struct_db3/
│   ├── pdb_snapshots/
│   └── extracted_pdbs/
└── runs/
```



## Example commands

```bash
mkdir -p inputs/resfiles protocols scripts outputs/ddg_db3 outputs/score_sc outputs/struct_db3 outputs/pdb_snapshots runs

sbatch --array=1-2114 scripts/run_flexddg_array.slurm inputs/pdbs/MUT_EPI_ISL_18133029_23_4muts_cngsug.pdb
python scripts/summarize_flexddg.py
bash scripts/extract_structures.sh
python scripts/rosetta_to_amber.py --ref reference.pdb
```
