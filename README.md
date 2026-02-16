# Cluster 29 Proliferation Analysis – scRNA-seq


This repository contains the analysis workflow used to identify and validate a proliferative cell population (Cluster 29) in a human single-cell RNA-seq dataset.

The analysis focuses on verifying whether Cluster 29 represents a proliferative state based on a curated cell-cycle gene signature and statistically significant differential expression results.

---

## Biological Context

The gene signature used in this analysis includes canonical proliferation and mitotic markers such as:

- Cell cycle regulators: CCNB2, CCNA2, PLK1, FOXM1, CDKN3
- Chromosome/kinetochore proteins: CENPH, NUF2, NCAPG
- Mitotic machinery: KIF2C, CEP55, UBE2C
- Proliferation markers: MKI67, BIRC5, TOP2A, TPX2

The objective is to determine whether Cluster 29 significantly expresses this proliferation program.


## Methods

1. Differential expression results were filtered using:
   - Adjusted p-value (p_val_adj < 0.05)
   - Genes belonging to the predefined proliferation signature

2. A proliferation score was computed per cell as the mean expression of significant signature genes.

3. Visualization:
   - UMAP 2D projection
   - UMAP 3D interactive plot
   - Violin plots by cluster
   - Cluster-level summary statistics
