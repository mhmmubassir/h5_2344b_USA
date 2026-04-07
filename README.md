# Phenotypic flexibility, rather than optimization, governs pandemic potential

Repository for the analyses, processed data products, figure-generation workflows, and supplementary materials associated with the manuscript:

**Mubassir et al. (2026)**
*Phenotypic flexibility, rather than optimization, governs pandemic potential*

---

## Summary

This repository contains the computational materials used to support a study of phenotypic diversification during the recent North American expansion of highly pathogenic avian influenza A H5N1 clade 2.3.4.4b. The project examines how phenotypic breadth emerges, persists, and is subsequently constrained across viral lineages, with a particular focus on hemagglutinin (HA), the major viral surface glycoprotein that mediates receptor binding, host range, and antigenically relevant structural variation.

Rather than treating adaptation as a simple process of progressive optimization toward a single best phenotype, this study evaluates whether pandemic-relevant emergence is better understood through the lens of **phenotypic flexibility**: the capacity of a lineage to explore and transiently maintain a broader region of phenotype space before ecological or host-associated filtering narrows that diversity. To address this question, the analyses integrate sequence-based, structure-based, and function-linked measurements across a large HA dataset spanning recent North American H5N1 circulation.

The repository therefore serves two purposes. First, it documents the analytical framework used in the manuscript, from sequence curation through structural modeling, phenotype extraction, and downstream statistical analyses. Second, it provides the processed outputs and figure-building components needed to reproduce the principal results and supplementary materials.

---

## Study scope

The analyses in this repository focus on the recent North American expansion of H5N1 clade 2.3.4.4b and evaluate phenotypic diversity across HA variants using a multi-layered framework that includes:

* **Sequence-derived phenotypes**, including biochemical and compositional descriptors
* **Structure-derived phenotypes**, including global and receptor-binding-site geometry metrics
* **Stability-oriented measurements**, derived from Rosetta-based calculations
* **Function-linked estimates**, including receptor-binding and glycan interaction analyses
* **Time-resolved diversity summaries**, used to evaluate expansion, pruning, and lineage-specific phenotypic breadth

Together, these analyses are used to test whether lineages associated with host shifts and broader cross-species potential are characterized less by convergence on a single optimal state and more by the transient maintenance of broader, more flexible phenotype repertoires.

---

## Repository contents

This repository contains the materials used to generate the main analyses presented in the manuscript, including:

* processed phenotype tables
* figure source data
* scripts for analysis and visualization
* figure panels and supplementary figures
* supplementary tables and associated outputs
* pipeline components used for structural and phenotypic characterization

Where appropriate, intermediate files have been reduced to processed or summarized outputs to improve portability and repository usability.

---

## Repository structure

```text
final_science_materials/
├── data/
├── main_figures/
├── pipelines/
├── supplementary_figures/
└── supplementary_tables/
```

---

## Analytical framework

The study integrates several complementary layers of analysis.

### 1. Sequence curation and variant definition

HA sequences were curated, aligned, and filtered to generate a high-confidence dataset suitable for temporal, structural, and phenotypic analysis. Unique amino acid variants were identified to reduce sampling redundancy while preserving biologically meaningful diversity across circulating lineages.

### 2. Structural modeling and geometry-based phenotyping

Representative HA variants were structurally modeled and analyzed to quantify both global structural divergence and local variation in the receptor-binding site. Particular emphasis was placed on loop geometry and subtle conformational shifts that may influence receptor engagement without requiring large-scale structural rearrangement.

### 3. Stability and energetic analyses

Rosetta-based calculations were used to estimate mutation-associated energetic effects and to benchmark computational predictions against experimentally informed measurements where available. These analyses were used not as absolute surrogates of function, but as one component of a broader phenotype framework.

### 4. Dynamics and glycan interaction analyses

Selected variants were further evaluated through molecular dynamics and post-simulation analysis workflows to examine receptor-binding-site flexibility, glycan-associated behavior, and linkage-dependent interaction patterns.

### 5. Diversity and lineage-level interpretation

Phenotypes were summarized through time and across lineages to test whether diversification precedes narrowing, and whether lineages associated with expanded host opportunity occupy broader phenotype space than lineages that remain closer to ancestral or avian-associated baselines.

---

## Reproducibility and workflow design

This repository is organized around reproducible, analysis-oriented outputs rather than raw high-volume runtime products. In most cases, the included scripts operate on processed or structured inputs and are intended to regenerate summary tables, visualizations, or final analytical products.

Because some upstream workflows depend on large external resources, high-performance computing environments, licensed software, or database-restricted sequence inputs, not every raw intermediate is distributed directly in the repository. Where that applies, the repository prioritizes:

* figure source data
* processed phenotype matrices
* reusable analysis scripts
* organized workflow components
* documentation sufficient to understand how major outputs were generated

Users interested in rerunning the complete upstream workflows should consult the pipeline-specific directories and local documentation for environment requirements and expected inputs.

---

## Software and computational environment

Different parts of the project rely on different computational environments. Depending on the workflow, analyses may require combinations of:

* Python
* shell scripting
* R or plotting libraries where applicable
* Rosetta
* AlphaFold-related modeling workflows
* Amber / cpptraj / MM/GBSA toolchains
* standard scientific Python libraries for statistics and figure generation

Pipeline-specific scripts may assume execution on a Linux-based high-performance computing environment with scheduler support.

---

## Data availability

This repository contains processed data products and manuscript-associated materials required to interpret and reproduce the reported analyses at the level of figures, tables, and summarized computational outputs.

Some upstream inputs, particularly those derived from external sequence repositories or large-scale raw structural simulations, may be subject to source-specific access conditions, file-size limitations, or separate archival practices. In such cases, this repository provides the processed derivatives necessary for the manuscript analyses while preserving a manageable and transparent project structure.

---

## Intended use

This repository is intended for:

* reproduction of manuscript figures and summary analyses
* inspection of phenotype-processing workflows
* reuse of downstream plotting and benchmarking utilities
* extension of the framework to related influenza HA datasets
* methodological reference for phenotype-centered viral emergence studies


## Public interactive tree

The public interactive tree is available at:

```text
https://nextstrain.org/community/mhmmubassir/h5_2344b_USA@main/full
```
---

## Citation

If you use this repository or any of its processed outputs, please cite:

**Mubassir et al. (2026)**
*Phenotypic flexibility, rather than optimization, governs pandemic potential*


---

## Contact

For questions regarding repository organization, computational workflows, or manuscript-associated materials, please contact the corresponding author.

---

## Notes

This repository reflects the analysis framework used for the manuscript and may continue to receive organizational improvements, documentation updates, and cleaned workflow components as the project is finalized for broader public release.
