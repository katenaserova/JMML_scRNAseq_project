# JMML scRNA-seq Project

# Single-cell transcriptomic analysis of JMML transformation toward AML

Master's thesis project

Author: Ekaterina Serova
ITMO University
Bioinformatics and Systems Biology

## Overview

This repository contains the code used for the Master's thesis project:


The project investigates cellular and molecular mechanisms associated with transformation of juvenile myelomonocytic leukemia (JMML) into acute myeloid leukemia (AML) using single-cell RNA sequencing (scRNA-seq) data.

The analysis includes:

* preprocessing and integration of scRNA-seq datasets
* cell type annotation
* identification of malignant populations using CopyKAT
* differential expression analysis
* GO and KEGG enrichment analysis
* pathway activity inference using PROGENy
* cell–cell communication analysis using CellChat

---

## Dataset

Due to ethical restrictions, raw sequencing data are not publicly available.

The study includes:

* 3 JMML-to-AML transformation bone marrow samples
* 2 reference JMML samples

---

## Software

Analyses were performed using:

* R (version 4.5.2)
* Seurat v5
* Harmony
* CopyKAT
* CellChat
* PROGENy
* clusterProfiler
* decoupleR

---

## Repository structure

### integration.Rmd

Integration of scRNA-seq datasets using Seurat v5 and Harmony.

Main steps:

* quality control
* normalization
* variable feature selection
* dimensionality reduction
* Harmony integration
* clustering
* UMAP visualization

---

### patient 1.Rmd

Patient-specific analysis for transforming JMML sample 1.

Includes:

* annotation review
* visualization
* exploratory analysis

---

### patient 2.Rmd

Patient-specific analysis for transforming JMML sample 2.

---

### patient 3.Rmd

Patient-specific analysis for transforming JMML sample 3.

---

### aneuploidy analysis.R

CopyKAT analysis.

Includes:

* inference of copy number alterations
* classification of cells as:

  * aneuploid
  * diploid
  * not defined

---

### DE and enrichment.R

Differential expression and enrichment analyses.

Includes:

* differential expression testing
* GO enrichment analysis
* KEGG enrichment analysis
* visualization of enrichment results

---

### cellchat_progeny.Rmd

CellChat and PROGENy analyses.

Includes:

* ligand–receptor interaction analysis
* communication network reconstruction
* pathway activity inference
* heatmap generation
* biological interpretation

---

## Reproducibility

To reproduce the analysis:

1. Install required R packages.
2. Place the processed Seurat objects in the working directory.
3. Run scripts in the following order:

```text
integration.Rmd
aneuploidy analysis.R
DE and enrichment.R
cellchat_progeny.Rmd
```

Patient-specific notebooks can be executed independently.


