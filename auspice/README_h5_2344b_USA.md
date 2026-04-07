# h5_2344b_USA

Interactive Nextstrain/Auspice phylogeny for U.S. H5N1 clade 2.3.4.4b HA sequences with sequence- and structure-informed phenotype overlays.

## Overview

This repository hosts an interactive dated phylogeny built from a BEAST-derived timetree and a 13,000-strain phenotype metadata table. The dataset is configured for visualization in Auspice / Nextstrain, allowing users to:

- color the tree by selected phenotypes
- inspect per-strain metadata
- view all variants across the dataset
- highlight the manuscript-updated top 10 variant classes

This build is intended for exploratory visualization and presentation.

## What this build shows

### Tree
- Time-scaled HA phylogeny for U.S. H5N1 clade 2.3.4.4b strains

### Variant views
- **All variants**: full variant labeling across the dataset
- **Top 10 variants**: manuscript-updated top 10 numbering shown with counts in parentheses

### Phenotype views
The tree can be colored by selected phenotype traits, including sequence-derived and structure-derived metrics.

To improve visualization, continuous phenotype color scales are clipped to the **2nd–98th percentile range**. This keeps extreme outliers from dominating the color range while preserving the underlying tree and metadata.

## Current visualization choices

### Variant labeling
- Global variant labels use the updated manuscript numbering where applicable
- Top 10 labels are shown in the updated numbering scheme with counts in parentheses

### Continuous color scaling
- Continuous phenotype fields are visualized using **2–98% clipped values**
- This affects color scaling only, not tree topology or strain membership

### Excluded traits
The following traits were intentionally removed from the current visualization because they were judged to add more noise than signal in this context:

- `flexddg_bind_23`
- `flexddg_bind_26`


## Public interactive tree

The public interactive tree is available at:

```text
https://nextstrain.org/community/mhmmubassir/h5_2344b_USA@main/full
```
## Maintainer

M H M Mubassir

