# Rosetta benchmarking with wet-lab validation

This module reduces the benchmarking branch to four scripts:

1. `1_group_ddg_outputs.py`
2. `2_summarize_rosetta_ddg.py`
3. `3_prepare_wetlab_benchmark_data.py`
4. `4_benchmark_rosetta_vs_wetlab.py`

## Recommended folder layout

```text
rosetta_benchmark/
├── scripts/
│   ├── 1_group_ddg_outputs.py
│   ├── 2_summarize_rosetta_ddg.py
│   ├── 3_prepare_wetlab_benchmark_data.py
│   └── 4_benchmark_rosetta_vs_wetlab.py
├── collected_ddg_outputs/
├── organized_ddg_groups/
├── bound_vs_unbound_heatmaps/
└── statistical_correlation/
    ├── experimental_ddg.csv
    ├── *.csv
    └── results/
```

## What each step does

### 1. Group raw output folders

```bash
python scripts/1_group_ddg_outputs.py \
  --src collected_ddg_outputs \
  --dst organized_ddg_groups
```

This creates grouped Rosetta/Flex ddG folders such as:

- `9dip_bn16_bb1/`
- `9dio_bn16_bb1_ori/`
- `9dio_bn16_bb1_rvmut/`
- `flexddg_all/`

### 2. Summarize Rosetta `.ddg` files

```bash
python scripts/2_summarize_rosetta_ddg.py \
  --groups organized_ddg_groups \
  --heatmaps bound_vs_unbound_heatmaps
```

This step:

- extracts avg and min ΔΔG from Rosetta `.ddg` files
- writes summary CSV/TXT tables
- creates heatmaps
- computes final bound-vs-HA-only ΔΔG tables

### 3. Prepare the benchmarking folder

```bash
python scripts/3_prepare_wetlab_benchmark_data.py \
  --search-root . \
  --outdir statistical_correlation
```

This step:

- generates `experimental_ddg.csv` from the embedded KD values
- collects final Rosetta benchmarking CSVs
- collects Flex ddG summary CSVs

### 4. Benchmark against wet-lab ΔΔG

```bash
python scripts/4_benchmark_rosetta_vs_wetlab.py \
  --indir statistical_correlation
```

Outputs go to:

- `statistical_correlation/results/`

including:

- `stats_summary.tsv`
- `stats_overall.csv`
- `all_predictions_merged.csv`
- overview scatter plots
- correlation bar plots
- overall performance plot

## Python packages

Required packages:

- `pandas`
- `numpy`
- `matplotlib`
- `scipy`
- `seaborn`

