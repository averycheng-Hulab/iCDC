# iCDC

---

# Single-Cell RNA-seq Analysis

This directory contains all scripts used for the single-cell RNA-seq workflows in this project, including whole-heart analysis, immune-cell profiling, subclustering, ECM analysis, and enrichment modules.

---

# 📁 Directory Structure

```

scRNA/
├── README.md                     # This file
│
├── utils/
│   └── utils_scRNA.R                   # Shared functions for all scRNA workflows
│
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
│   ├── cd45_b_subcluster.R             # B-cell subclustering
│   └── cd4_T_TCR_module.R              # TCR clonotype + CD4 subset integration module
│
├── ECM_genelist/
│   └── ECMregulators.csv           # Gene list for ECM regulatory module score calculation
│
└── Enrichment/
├── enrich_scRNA.R                  # Unified IR vs Sham/DC enrichment (any celltype)
└── enrich_scRNA_cd4.R              # CD4-specific MI vs iCDC enrichment

````

---

# **1. Utilities**

### **utils_scRNA.R**

This file provides a unified set of helper functions used across all scRNA pipelines:

* SCTransform → PCA → Harmony integration wrapper  
* Clustering workflow (Neighbors, FindClusters, UMAP, t-SNE)  
* QC pipeline with DoubletFinder  
* Cell-cycle scoring wrapper  
* Filtering helpers (mitochondrial / ribosomal genes)  
* Marker detection wrapper  
* Annotation mapping utilities  
* Unified ggplot theme + save functions  
* Batch-agnostic processing for all Seurat objects  

All analysis scripts use:

```r
source("utils_scRNA.R")
````

---

# **2. Main Pipelines**

## **whole_heart_main.R**

Full whole-heart workflow:

* Read all 10× matrices
* QC + DoubletFinder
* Cell-cycle scoring
* SCTransform → PCA → Harmony integration
* Clustering + UMAP/tSNE
* Save processed object for downstream modules

---

## **cd45_main.R**

CD45+ immune compartment workflow:

* Extract immune cells
* Re-integration
* Clustering + annotation
* Save object for subclustering modules

---

## **cd4_T_main.R**

CD4 T cell workflow:

* Subset CD4 T cells
* Integration and clustering
* Activation, cytokine, regulatory module scoring
* Saves object for CD4-TCR module & CD4 enrichment module

---

# **3. Subcluster Pipelines**

Located under `scRNA/Subcluster/`.

---

## **Major fibroblast modules**

### **whole_heart_fibroblast.R**

* Subset fibroblasts
* Re-integrate & re-cluster
* Annotate fibroblast subtypes (`FibType`, `FibType3`)

### **whole_heart_fibroblast_ecm.R**

**ECM regulator-only analysis**

Computes:

* **ECM_Regulator_Score** 

Outputs:

* UMAP visualization
* FibType violin plots
* Treatment violin plots
* Saves updated Seurat object

---

## **Immune subclustering modules**

All follow the re-integration → re-cluster → refine → annotate pattern.

### **cd45_MacMonoDc_subcluster.R**

Monocyte/Macrophage/DC subclustering with multi-step refinement.

### **cd45_tnk_subcluster.R**

T/NK subclustering and supervised annotation.

### **cd45_neutrophil_subcluster.R**

Two-round integration, removal of unwanted clusters, annotation.

### **cd45_b_subcluster.R**

B-cell subclustering:

* Re-integrate
* Annotate into B-Naive / B-Act / Memory / Plasma

---

## **cd4_T_TCR_module.R** 

Integration of **TCR clonotype** information into CD4 T-cell subclusters.

This module:

* Reads CD4+ T-cell Seurat object
* Imports TCR clonotype tables (filtered contig annotations)
* Merges clonotype metadata into Seurat metadata
* Identifies clonal expansion patterns within CD4 subsets
* Produces:

  * Clonal expansion UMAP
  * Clonotype frequency tables
  * TCR gene usage analysis (TRAV/TRBV patterns)
* Saves enriched CD4 object with TCR info

This script connects the CD4 phenotypic clusters with their clonal architecture.

---

# **4. Enrichment Modules**

Located in `scRNA/Enrichment/`.

---

## **enrich_scRNA.R**

General enrichment module for *any* scRNA cell type:

* DEG: IR vs Sham
* DEG: IR vs DC
* Intersect shared up/down DEGs
* GO enrichment
* KEGG enrichment
* Bubble plots
* Excel export

Usage example:

```r
run_enrich_scRNA(obj, celltype = "F-SL", outdir = "FSL_enrich")
```

---

## **enrich_scRNA_cd4.R**

CD4-specific enrichment:

* DEG: MI vs iCDC
* GO/KEGG enrichment
* Bubble plots
* Organizes output into CD4-specific folder structure

---


