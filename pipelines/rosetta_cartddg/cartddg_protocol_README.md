# Cartesian ddG protocol layout

Suggested repo layout:

```text
cartddg_protocol/
├── protocols/
│   └── cart_relax.script
├── inputs/
│   ├── pdbs/
│   └── mutfiles/
├── scripts/
│   ├── 1_cart_relax.slurm
│   ├── 2_sort_relax_scores.py
│   ├── 3_select_best_relaxed_model.py
│   └── 4_cartddg_array.slurm
├── outputs/
│   ├── cart_relax/
│   │   ├── models/
│   │   ├── scores/
│   │   └── best_models/
│   └── cartddg/
│       ├── ddg/
│       ├── scores/
│       ├── wt_pdbs/
│       └── mut_pdbs/
├── runs/
└── logs/
```

## Workflow

1. Run Cartesian relax on the starting glycoprotein PDB.
2. Sort the relax score file.
3. Select the best relaxed model.
4. Run Cartesian ddG as a Slurm array over all `.mut` files.

## Example commands

```bash
sbatch scripts/1_cart_relax.slurm inputs/pdbs/your_input.pdb
python scripts/2_sort_relax_scores.py outputs/cart_relax/scores/your_input_relax.sc
python scripts/3_select_best_relaxed_model.py \
    outputs/cart_relax/scores/your_input_relax.sc \
    outputs/cart_relax/models/your_input \
    -o outputs/cart_relax/best_models/your_input_relaxed_best.pdb
sbatch --array=1-N scripts/4_cartddg_array.slurm \
    outputs/cart_relax/best_models/your_input_relaxed_best.pdb
```
