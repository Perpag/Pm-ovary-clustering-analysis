
**Patiria miniata adult ovary scripts**

This repository contains scripts used for the manuscript "Molecular evidence for early deuterostome origins of ovarian cell types and neuroendocrine control of reproduction " by Periklis Paganos, Beverly Naigles, Carsten Wolff, Danila Voronov, S. Zachary Swartz.


# System Requirements

## Operating System

-   macOS 12+, Linux (Ubuntu 20.04+), Windows 10+
-   Tested on macOS Sonoma, macOS Sequoia, macOS 26 Tahoe

## Software Dependencies (with tested versions)

### R Environment

-   **R ≥ 4.4.1**
-   **Required R packages:**
    -   Seurat ≥ 5.3.0
    -   dplyr ≥ 1.1
    -   ggplot2 ≥ 3.4
    -   Matrix ≥ 1.5
    -   SAMap (Python package), via reticulate or precomputed matrices
    -   reticulate ≥ 1.40.0
    -   plotly ≥ 4.10
    -   pheatmap ≥ 1.0
    -   tidyverse ≥ 2.0.0
    -   cowplot ≥ 1.1.3
    -   monocle3 ≥ 1.3.7
    -   SeuratWrappers ≥ 0.3.5
    -   patchwork ≥ 1.3.0
    -   EnhancedVolcano ≥ 1.22.0
    -   reshape2 ≥ 1.4.4
    -   ggalluvial ≥ 0.12.5
    -   stringr ≥ 1.5.1
    -   ggnewscale ≥ 0.5.2
    -   ggrepel ≥ 0.0.6
    -   SeuratData ≥ 0.2.2.9001
    -   SeuratDisk ≥ 0.0.0.9021
    -   topGO ≥ 2.56.0
    -   clusterProfiler ≥ 4.12.6
    -   networkD3 ≥ 0.4

### Python

-   Python ≥ 3 (required for SAMap)
-   SAMap ≥ 1.0.15

## Hardware Requirements

A typical laptop or workstation with sufficient RAM to support in-memory
operations.

------------------------------------------------------------------------

# Installation Guide

## 1. Install R

Download and install R from CRAN.

## 2. Install Required R Packages

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages(c(
  "pheatmap", "tidyverse", "cowplot", "dplyr", "patchwork",
  "reshape2", "ggalluvial", "ggplot2", "stringr", "ggnewscale",
  "ggrepel", "networkD3", "plotly"
))

install.packages(c("Matrix", "reticulate"))
install.packages("Seurat")
install.packages("SeuratData")
install.packages("SeuratDisk")
install.packages("EnhancedVolcano")

BiocManager::install("monocle3")
BiocManager::install(c("clusterProfiler", "topGO"))
```

## 3. Install Python + SAMap

Follow installation instructions from the SAMap repository: -
https://github.com/atarashansky/SAMap

## 4. Download Datasets

Download the datasets from NCBI GEO (accession: **GSE292963**).

### Typical Install Time

-   R packages: **2--4 minutes**
-   Python + SAMap: **2--4 minutes**

------------------------------------------------------------------------

# Demo

## Instructions to Run the Analysis

1.  Visit the dataset GEO page:
    -   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292963
2.  Place the files into your working directory.
3.  Run the script:
    -   https://github.com/Perpag/Pm-ovary-clustering-analysis/blob/main/Clustering_analysis_pm_ovary.R

## Expected Output

-   Clustering of *P. miniata* ovary single-cell data
-   UMAP reconstruction

## Expected Run Time

-   Full dataset: **\~20--60 minutes** (depending on system and dataset
    size)

------------------------------------------------------------------------

# Full Workflow Instructions

## Clustering Analysis (R)

The `Clustering_analysis_pm_ovary.rmd` file contains the full workflow
for analyzing *P. miniata* ovary single-cell data. It uses: - **Seurat /
SeuratWrappers** for preprocessing, normalization, integration,
clustering, PCA, UMAP. - **Visualization packages:** ggplot2, cowplot,
patchwork, pheatmap. - **Utility + wrangling:** reshape2, tidyverse,
tidyr, dplyr. - **Matrix operations:** matrixStats. - **Trajectory
analysis:** monocle3.

## Cross-Species Comparison with SAMap

The `pm_cross-species_comparison_with_SAMap.rmd` file explains R--Python
integration for cross-species single-cell atlas comparison.

### R Environment

-   Seurat, SeuratData, SeuratDisk
-   reshape2 for data reshaping
-   ggalluvial + reshape2 for sankey-style figures

### Python Environment

-   **samap**
-   Tools: get_mapping_scores, GenePairFinder, sankey_plot, chord_plot,
    CellTypeTriangles, ParalogSubstitutions, FunctionalEnrichment,
    convert_eggnog_to_homologs, GeneTriangles

## Gene Ontology Enrichment

The `pm_gene_ontology_enrichment.rmd` file describes GO enrichment
analysis and visualization of DE genes across germline clusters.

## Additional Dataset Analyses

The `Clustering_analysis_of_hs_fetal_ovary_and_mm_pituitary_gland.R`
script analyzes human fetal ovary and mouse pituitary datasets, using
the same Seurat-based workflow.

## Cell Ranger Mapping

The `pm_cellranger_non_forced.sh` script shows how to map sequencing
libraries to the reference genome.

------------------------------------------------------------------------

# Pipeline Overview

-   Preprocessing
-   Clustering
-   SAMap cross-species mapping
-   Visualization
-   Export of results
