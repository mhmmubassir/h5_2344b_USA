# Phenotype workflow

This organized set places the scripts in analytical order.

1. `1_calculate_sequence_phenotypes.py`
   Computes sequence-derived HA phenotype metrics and a reference-relative delta table.

2. `2_compute_phenotypic_diversity_allseq_bins5.py`
   Uses phenotype tables to quantify phenotypic diversity through time using 5-bin Shannon effective bins with bootstrap/rarefaction.

3. `3_summarize_lineage_A_B_phenotypes.py`
   Summarizes phenotype distributions by lineage group using median and IQR tables.

4. `4_plot_tree_tip_colors_v3_v10_anchors.py`
   Colors phylogenetic tree tips using precomputed V3/V10 anchor-based lineage-likeness scores.
