# iCDC


---



# Single-Cell RNA-seq Analysis (scRNA)

This directory contains all scripts used for the single-cell RNA-seq workflows in this project, including whole-heart analysis, immune-cell profiling, subclustering, ECM analysis, and enrichment modules.

---

# 📁 Directory Structure

```
scRNA/
├── README_scRNA.md                     # This file
│
├── utils/
│   └── utils_scRNA.R                   # Shared functions for all scRNA workflows
│
├── plots/                              # (Empty)
├── results/                            # (Empty)
│
├── Main/
│   ├── whole_heart_main.R              # Whole-heart scRNA integration workflow
│   ├── cd45_main.R                     # CD45+ immune compartment analysis
│   └── cd4_T_main.R                    # CD4+ T-cell analysis
│
├── Subcluster/
│   ├── whole_heart_fibroblast.R        # Fibroblast re-clustering
│   ├── whole_heart_fibroblast_ecm.R    # Fibroblast ECM regulator program analysis
│   ├── cd45_MacMonoDc_subcluster.R     # Monocyte / Macrophage / DC subclustering
│   ├── cd45_tnk_subcluster.R           # T/NK subclustering
│   ├── cd45_neutrophil_subcluster.R    # Neutrophil subclustering
│   └── cd45_b_subcluster.R             # B-cell subclustering
│
└── Enrichment/
    ├── enrich_scRNA.R                  # Unified IR vs Sham/DC enrichment (any celltype)
    └── enrich_scRNA_cd4.R              # CD4-specific MI vs iCDC enrichment
```

---

#  **1. Utilities**

### **utils_scRNA.R**

This file provides a unified set of helper functions used across all pipelines:

* SCTransform → PCA → Harmony integration wrapper
* Clustering (Neighbors, FindClusters, UMAP, t-SNE)
* Sample-wise QC with DoubletFinder
* Cell-cycle scoring wrapper
* Filtering helpers (mitochondrial / ribosomal / hemoglobin genes)
* Marker detection wrapper
* Annotation mapping utilities
* Plotting functions with a consistent theme
* Batch-agnostic processing for all major Seurat objects

All main and subcluster scripts source this file:

```r
source("utils_scRNA.R")
```

---

#  **2. Main Pipelines**

## **whole_heart_main.R**

Performs whole-heart scRNA-seq analysis, including:

* Reading all 10× matrices
* QC + DoubletFinder
* Cell-cycle scoring
* SCTransform → PCA → Harmony
* Clustering and visualization
* Saving processed whole-heart object

## **cd45_main.R**

Extracts CD45+ immune compartment and performs:

* Integration
* Clustering
* Annotation
* Stores object for downstream subclustering

## **cd4_T_main.R**

Focused CD4 T cell workflow:

* Subsetting
* Integration + clustering
* T cell activation/regulatory module scoring
* Output used by CD4 enrichment modules

---

#  **3. Subcluster Pipelines**

Located under `scRNA/Subcluster/`.

## **Major fibroblast modules**

### **whole_heart_fibroblast.R**

* Extracts fibroblast clusters
* Re-integrates and re-clusters
* Annotates fibroblast subtypes (`FibType`, `FibType3`)

### **whole_heart_fibroblast_ecm.R**

 **ECM regulator-only analysis (per your request)**
Computes **only**:

* **ECM_Regulator_Score**

Then visualizes:

* UMAP expression
* Violin plots by FibType
* Violin plots by Treatment
* Saves updated fibroblast object

No collagen / glycoprotein / proteoglycan / combined ECM modules are calculated.

---

## **Immune subclustering modules**

All follow the same pattern:

### **cd45_MacMonoDc_subcluster.R**

Re-clusters macrophage/monocyte/DC compartments, includes multi-step refinement.

### **cd45_tnk_subcluster.R**

Subclustering of T/NK modules.

### **cd45_neutrophil_subcluster.R**

Two-round integration → removal of unwanted clusters → supervised annotation.

### **cd45_b_subcluster.R**

B-cell subset analysis:

* Re-integration
* Supervised annotation into B-Naive / B-Act / Memory / Plasma
* Marker tables

---

#  **4. Enrichment Modules**

Located under `scRNA/Enrichment/`.

## **enrich_scRNA.R**

General enrichment module for any scRNA celltype:

* DEG: IR vs Sham
* DEG: IR vs DC
* Intersection of up/down DEGs
* GO enrichment
* KEGG enrichment (MMU)
* Bubble plots
* Excel export (all sheets in one file)

Usage:

```r
run_enrich_scRNA(obj, celltype = "F-SL", outdir = "FSL_enrich")
```

---

## **enrich_scRNA_cd4.R**

CD4 T-cell specific enrichment:

* DEG: MI vs iCDC
* GO / KEGG enrichment
* Bubble plot visualization
* CD4-specific output folder structure

---

#  **5. Data Availability and Privacy**

Because this repository corresponds to **unpublished work**, all folder locations intended to contain raw or processed data are left intentionally empty:

* `scRNA/plots/`
* `scRNA/results/`
* `bulk/input/`
* `bulk/output/`
* `.Rds` files required as input

Only scripts are provided.

> **To protect unpublished data, all data and result folders are intentionally left empty.
> Only code is provided to ensure reproducibility without exposing sensitive information.**

---

#  **6. Software Requirements**

* R ≥ 4.1
* Seurat ≥ 4.3
* harmony
* DoubletFinder
* clusterProfiler
* org.Mm.eg.db
* writexl
* ggplot2 / patchwork / dplyr / stringr

---

#  **7. Maintainer: Guo Cheng**


---
