# FoldX protocol

Suggested serial workflow:

1. `1_foldx_repair.slurm`
2. `2_foldx_buildmodel_ph_array.slurm`
3. `3_foldx_summarize_ddg.py`

## Suggested repo layout

```text
foldx_protocol/
├── inputs/
│   ├── pdbs/
│   └── foldx_mutfiles/
├── scripts/
│   ├── 1_foldx_repair.slurm
│   ├── 2_foldx_buildmodel_ph_array.slurm
│   └── 3_foldx_summarize_ddg.py
├── outputs/
│   ├── foldx_repair/
│   ├── foldx_buildmodel/
│   │   └── per_variant/
│   └── foldx_summary/
└── runs/
```

## Step 1: Repair the input structure

```bash
sbatch scripts/1_foldx_repair.slurm inputs/pdbs/MUT_EPI_ISL_18133029_4muts.pdb
```

Optional environment overrides:

```bash
ION=0.05 TEMP=298 REPAIR_PH=7.0 sbatch scripts/1_foldx_repair.slurm inputs/pdbs/model.pdb
```

Output:

- `outputs/foldx_repair/<input>_Repair.pdb`
- `outputs/foldx_repair/*.fxout`
- `outputs/foldx_repair/<input>_repair_metadata.tsv`

## Step 2: Run BuildModel across pH values

Use the repaired PDB from step 1.

```bash
sbatch --array=1-N scripts/2_foldx_buildmodel_ph_array.slurm outputs/foldx_repair/MUT_EPI_ISL_18133029_4muts_Repair.pdb
```

Optional environment overrides:

```bash
MUTDIR=inputs/foldx_mutfiles \
NRUNS=3 \
TEMP=298 \
ION=0.05 \
PH_LIST_CSV=4.5,5.0,5.3,5.5,5.7,6.0,6.5,7.0 \
sbatch --array=1-N scripts/2_foldx_buildmodel_ph_array.slurm outputs/foldx_repair/model_Repair.pdb
```

Per-variant outputs are organized as:

```text
outputs/foldx_buildmodel/per_variant/<variant>/T298/
├── fxout/
├── failures.tsv
├── individual_list.txt
└── run_metadata.tsv
```

## Step 3: Summarize FoldX outputs

```bash
python3 scripts/3_foldx_summarize_ddg.py \
  --root outputs/foldx_buildmodel/per_variant \
  --outdir outputs/foldx_summary
```

Optional manifest mapping mutation-file names to sequence IDs or tree taxa:

```bash
python3 scripts/3_foldx_summarize_ddg.py \
  --manifest inputs/foldx_mutfiles/foldx_mutfiles_manifest.tsv
```
