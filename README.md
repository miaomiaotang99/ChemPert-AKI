# AKI Single-Cell Virtual Perturbation Framework

## Overview

This repository provides a reproducible computational framework for investigating chemical perturbation responses in acute kidney injury (AKI) by integrating:

- bulk RNA-seq validation
- single-cell RNA-seq analysis
- immune infiltration analysis
- cell-cell communication analysis
- network-based regulatory analysis
- virtual gene perturbation
- molecular representation learning
- deep-learning-based perturbation prediction

The overall goal is to establish a computational pipeline for predicting cellular transcriptional responses induced by chemical exposure (e.g., DEHP/MEHP) and identifying molecular mechanisms associated with AKI progression.

---

# Repository Architecture

```
AKI-DEHP-scRNAseq/
│
├── bulk_analysis/
│   ├── GSE139061_normalization.R
│   ├── CIBERSORT_analysis.R
│   ├── immune_violin_plot.R
│   └── gene_immune_correlation_lollipop.R
│
├── single_cell_analysis/
│   ├── 01_Create_Seurat_Object.R
│   ├── 02_Normalization_BatchCorrection_DimensionReduction.R
│   ├── 03_SingleCell_Annotation.R
│   ├── 04_Enrichment_Analysis.R
│   ├── 05_Volcano_Plot.R
│   ├── 06_Cell_Interaction.R
│   ├── 07_CellChat_Communication.R
│   ├── 08_Virtual_Knockout_Analysis.R
│   └── 09_hdWGCNA_Analysis.R
│
├── perturbation_model/
│   ├── model/
│   └── fusion_cpa/
│
└── README.md
```

---

# Single-cell RNA-seq Workflow

## 1. Seurat Object Construction

Input:

- 10X Genomics count matrix

Processing:

- Create Seurat object
- Quality control
- Mitochondrial/ribosomal filtering
- Save processed object

Output:

```
seurat_raw_object.rds
```

---

## 2. Normalization and Batch Correction

Methods:

- Log normalization
- Highly variable gene selection
- PCA
- Harmony batch correction
- UMAP/tSNE visualization
- Leiden/Louvain clustering

Output:

```
seurat_processed_object.rds
```

---

## 3. Cell Type Annotation

Cell populations are annotated using:

- canonical marker genes
- DotPlot visualization
- UMAP/tSNE visualization

Major kidney cell types:

- Proximal tubule cells
- TAL cells
- DCT cells
- Collecting duct cells
- Podocytes
- Endothelial cells
- Immune cells
- Interstitial cells

---

# Functional Analysis

## Enrichment Analysis

Implemented:

- GO enrichment
- Hallmark GSEA
- KEGG enrichment
- Reactome enrichment

Output:

- pathway enrichment tables
- GSEA curves


---

# Immune Landscape Analysis

Bulk RNA-seq immune profiling:

Algorithm:

- CIBERSORT

Outputs:

- immune cell fractions
- Control vs AKI immune differences
- violin plots
- immune correlation analysis


---

# Cell Communication Analysis

## CellChat

Analysis includes:

- ligand-receptor interaction inference
- signaling pathway activity
- communication network visualization


---

# Network Analysis

## hdWGCNA

Functions:

- metacell construction
- co-expression network construction
- module identification
- hub gene detection
- module scoring


---

# Virtual Gene Perturbation

## scTenifoldKnk

Purpose:

Simulate:

- single-gene knockout
- combinatorial knockout

Examples:

```
STAT3
TP53
AKT1
```

Outputs:

- differential regulatory genes
- network perturbation effects


---

# Molecular Perturbation Prediction Model

## Baseline

This repository keeps chemCPA as the original baseline model.

chemCPA predicts:

> Cellular responses to unseen chemical perturbations at single-cell resolution.

The original framework uses molecular embeddings and perturbation-response modeling.

---

# chemCPA × Morgan–D-MPNN Fusion Model

## Overview

The fusion framework extends chemCPA by integrating:

- Morgan fingerprint representation
- directed message passing neural network (D-MPNN)
- cross-attention fusion
- dose encoding
- cell-state embedding
- Gaussian decoder


Architecture:

```
SMILES
 |
 +----------------+
 |                |
Morgan FP      D-MPNN
 |                |
 +----Fusion------+
        |
 Chemical latent z
        |
 Dose + Cell type + Basal expression
        |
 Gaussian Decoder
        |
 Predicted perturbation expression
```

---

# Fusion Modes

Four molecular representation strategies:

| Mode            | Description                           |
| --------------- | ------------------------------------- |
| morgan          | Morgan fingerprint only               |
| graph           | Graph neural network only             |
| concat          | Morgan + graph concatenation          |
| cross_attention | Morgan + graph cross-attention fusion |


---

# Model Input

Training CSV:

| Column         | Description                |
| -------------- | -------------------------- |
| smiles         | canonical molecular SMILES |
| dose           | perturbation concentration |
| cell_type      | cellular identity          |
| control_path   | control expression `.npy`  |
| perturbed_path | treated expression `.npy`  |


Expression format:

```
[G genes]
```

All samples must share:

- identical gene order
- identical normalization strategy


---

# Training

Example:

```bash
python -m fusion_cpa.train \
--csv data/train.csv \
--mode cross_attention \
--epochs 100 \
--batch-size 128
```

Prediction:

```bash
python -m fusion_cpa.predict \
--checkpoint model.pt \
--smiles "CCO"
```

---

# Benchmark Evaluation

Recommended metrics:

- Pearson correlation coefficient
- R²
- MSE/RMSE
- Energy distance
- Maximum mean discrepancy (MMD)

Evaluation settings:

- seen compounds
- unseen compounds
- unseen cell types
- dose extrapolation


---

# Data Availability

The repository does not include raw sequencing datasets.

Required external resources:

- GEO datasets
- PubChem chemical information
- CTD disease associations
- CIBERSORT signature matrix



