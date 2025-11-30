# iCDC
# Utilities for scRNA-seq Analysis
This folder contains `utils_scRNA.R`, a shared function library used across
all single-cell RNA-seq analysis modules, including:

- Main workflows (`whole_heart_main.R`, `cd45_main.R`, `cd4_T_main.R`)
- Subclustering modules (fibroblast, Mac/Mono/DC, TNK, Neutrophil, B-cell)
- Enrichment modules
- TCR module (indirectly, for Seurat object consistency)

The utilities ensure unified preprocessing, batch correction, clustering,
marker detection, and mapping logic across all scRNA pipelines.

---

## 📌 Contents of utils_scRNA.R

### 1. **Preprocessing Wrappers**
Reusable wrappers for standard single-cell preprocessing:

- `run_sct_pca_harmony()`
  - SCTransform (with optional regression)
  - RunPCA
  - Harmony batch correction
  - Returns updated Seurat object

- `run_clustering()`
  - FindNeighbors / FindClusters
  - UMAP + TSNE embeddings
  - Configurable dims + resolution

These wrappers ensure that all main and subcluster pipelines use the same
integration and clustering logic.

---

## 2. **Marker & DEG Utilities**
Functions used by main pipelines and subclustering:

- `find_all_markers_wrapper()`
- `save_markers()`
- tidy marker output formatter  
- auto-export to CSV

---

## 3. **Annotation Tools**
General annotation helper functions:

- `annotate_by_cluster()`  
  Apply a named vector mapping (cluster → label) to a Seurat object.

- `merge_annotation_back()`  
  Merge refined annotations (from a subcluster object) back into the
  parent object.

These tools are used extensively in:
- Fibroblast refinement (FibType1 / FibType2 / FibType)
- Mac/Mono/DC (MacType1 → refined → final)
- TNK / Neutrophil / B-cell supervised mappings

---

## 4. **Subclustering Support (ALL Modules)**

This section provides reusable helpers for all downstream subclustering modules.

### 🟦 **Fibroblast Subclustering**
Used in:
- `whole_heart_fibroblast.R`
- `whole_heart_fibroblast_ecm.R`

Functions support:
- extraction of fibroblast clusters  
- re-normalization and Harmony batch correction  
- ECM score generation  
- fibroblast type annotation merging  

### 🟧 **Macrophage / Monocyte / DC**
Used in `cd45_MacMonoDc_subcluster.R`:

- selection of immune-myeloid clusters  
- global + fine SCTransform reprocessing  
- supervised annotation (MacType1 / MacType2 / refined)  
- merge refined MacType into the CD45 main object  

### 🟩 **TNK subclustering**
Used in `cd45_tnk_subcluster.R`:

- selection of T/NK clusters  
- Harmony correction & re-clustering  
- TNKType1 annotation  

### 🟨 **Neutrophil subclustering**
Used in `cd45_neutrophil_subcluster.R`:

- extraction of neutrophils  
- clustering and annotation into NeutType1  

### 🟧 **B-cell subclustering**
Used in `cd45_b_subcluster.R`:

- clustering  
- supervised BType1 annotation  

All modules use the same SCTransform → PCA → Harmony → clustering utilities.

---

## 5. **Plotting Utilities**
Simple ggplot-based functions used across the project:

- common theme (`plot_theme_common()`)
- save function (`save_plot()`)
- bubble plot helpers (if needed by enrichment)

---

## 6. **Why This Matters**
Without `utils_scRNA.R`, each subpipeline would independently duplicate the
same methods—making the analysis inconsistent and harder to maintain.

With utilities:

- All batch correction is identical  
- All clustering parameters are standardized  
- All annotation mapping uses the same API  
- Subcluster modules become shorter and safer  
- Cross-module reproducibility is 100% guaranteed  

---

## Location
`scRNA/utils/utils_scRNA.R`

This file is sourced by all main and subcluster scripts:
```r
source("scRNA/utils/utils_scRNA.R")
