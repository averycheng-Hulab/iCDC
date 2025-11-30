###############################################################################
# Whole-heart scRNA-seq pipeline (10X, Seurat + DoubletFinder + Harmony)
#
# - Input: 12 filtered_feature_bc_matrix folders (DC / IR / Sham / Vector)
# - Steps (per sample): QC, filtering, DoubletFinder (singlets only)
# - Steps (merged): merge → remove mito/ribo → cell cycle scoring
#                   → SCTransform (regress CC.Difference)
#                   → PCA → Harmony → clustering → UMAP/TSNE
# - Output:
#     r_objects/seuobj_final.Rds
#     markers_by_cluster.csv
#     markers_by_CellType1.csv
###############################################################################

## ---------------------------------------------------------------------------
## 0. Load packages & utils
## ---------------------------------------------------------------------------

required_pkgs <- c(
  "Seurat","DoubletFinder","harmony","ggplot2","patchwork",
  "dplyr","stringr","clustree","scater","SCP","writexl"
)
missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse=", "))
})

library(Seurat)
library(DoubletFinder)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(clustree)
library(scater)
library(SCP)
library(writexl)

source("utils_scRNA.R")

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)

dir.create("QC", showWarnings = FALSE)
dir.create("output_plots", showWarnings = FALSE)
dir.create("r_objects", showWarnings = FALSE)


## ---------------------------------------------------------------------------
## 1. Sample names & paths
## ---------------------------------------------------------------------------

samples <- c(
  "DC_1","DC_2","DC_3",
  "IR_1","IR_2","IR_3",
  "Sham_1","Sham_2","Sham_3",
  "Vector_1","Vector_2","Vector_3"
)

treatment_map <- c(
  "DC_1"="DC","DC_2"="DC","DC_3"="DC",
  "IR_1"="IR","IR_2"="IR","IR_3"="IR",
  "Sham_1"="Sham","Sham_2"="Sham","Sham_3"="Sham",
  "Vector_1"="Vector","Vector_2"="Vector","Vector_3"="Vector"
)

base_dir <- "Matrix"
matrix_dirs <- file.path(base_dir, samples, "filtered_feature_bc_matrix")
stopifnot(length(samples) == length(matrix_dirs))


## ---------------------------------------------------------------------------
## 2. Clustering parameters
## ---------------------------------------------------------------------------

pc_dims     <- 1:40
cluster_res <- seq(0.2, 1.6, 0.2)


## ---------------------------------------------------------------------------
## 3. Per-sample QC + DoubletFinder
## ---------------------------------------------------------------------------

seuobj_list <- vector("list", length(samples))
names(seuobj_list) <- samples

for (i in seq_along(samples)) {
  seuobj_list[[i]] <- process_sample_scRNA(
    matrix_dir      = matrix_dirs[i],
    sample_name     = samples[i],
    treatment_label = treatment_map[[samples[i]]]
  )
}


## ---------------------------------------------------------------------------
## 4. Merge singlet-filtered samples
## ---------------------------------------------------------------------------

seuobj_merged <- seuobj_list[[1]]
if (length(seuobj_list) > 1) {
  for (i in 2:length(seuobj_list)) {
    seuobj_merged <- merge(seuobj_merged, y = seuobj_list[[i]])
  }
}

BasicData <- GetAssayData(seuobj_merged, assay="RNA", slot="counts")
meta_min  <- seuobj_merged@meta.data %>% select(orig.ident, Treatment)

seuobj_raw <- CreateSeuratObject(counts = BasicData, meta.data = meta_min)
seuobj_raw$orig.ident <- factor(seuobj_raw$orig.ident, levels = samples)
seuobj_raw$Treatment  <- factor(seuobj_raw$Treatment, levels = c("Sham","IR","DC","Vector"))


## ---------------------------------------------------------------------------
## 5. Remove mitochondrial & ribosomal genes
## ---------------------------------------------------------------------------

seuobj_qc <- filter_mito_ribo(seuobj_raw)


## ---------------------------------------------------------------------------
## 6. Cell cycle scoring
## ---------------------------------------------------------------------------

if (!exists("cc.genes.updated.2019")) {
  data("cc.genes.updated.2019", package="Seurat")
}

s_genes   <- stringr::str_to_title(cc.genes.updated.2019$s.genes)
g2m_genes <- stringr::str_to_title(cc.genes.updated.2019$g2m.genes)

seuobj_qc <- compute_cell_cycle_scores(
  seu       = seuobj_qc,
  s_genes   = s_genes,
  g2m_genes = g2m_genes
)


## ---------------------------------------------------------------------------
## 7. SCTransform → PCA → Harmony
## ---------------------------------------------------------------------------

seuobj_final <- run_sct_pca_harmony(
  seu             = seuobj_qc,
  vars_to_regress = "CC.Difference",
  npcs            = 50,
  dims_use        = pc_dims,
  harmony_group   = "orig.ident"
)


## ---------------------------------------------------------------------------
## 8. Clustering
## ---------------------------------------------------------------------------

for (res in cluster_res) {
  seuobj_final <- run_clustering(
    seu        = seuobj_final,
    dims_use   = pc_dims,
    resolution = res
  )
}

ctree <- clustree(seuobj_final@meta.data, prefix="SCT_snn_res.")
save_plot(ctree, "wholeheart_clustree.pdf", width=15, height=10)


## ---------------------------------------------------------------------------
## 9. Markers before annotation
## ---------------------------------------------------------------------------

markers_by_cluster <- find_all_markers_wrapper(
  seu   = seuobj_final,
  ident = "seurat_clusters"
)
write.csv(
  markers_by_cluster,
  file = "r_objects/markers_by_cluster.csv",
  row.names = FALSE
)


## ---------------------------------------------------------------------------
## 10. Annotation (CellType1)
## ---------------------------------------------------------------------------

cluster_to_celltype <- c(
  "0"="Fibroblast","1"="Fibroblast","2"="Fibroblast","3"="Fibroblast",
  "4"="Fibroblast","5"="Fibroblast","6"="Fibroblast","9"="Fibroblast",
  "12"="Fibroblast","13"="Fibroblast","17"="Fibroblast","25"="Fibroblast",
  "28"="Fibroblast","33"="Fibroblast","37"="Fibroblast","46"="Fibroblast",

  "45"="Mast Cell","31"="Erythrocyte",
  "27"="Macrophage","36"="Monocyte",
  "32"="cDC","44"="pDC",

  "10"="Endothelium","11"="Endothelium","18"="Endothelium",
  "19"="Endothelium","20"="Endothelium","24"="Endothelium","40"="Endothelium",

  "34"="Cycling Cell",

  "7"="SMC_Pericyte","22"="SMC_Pericyte","35"="SMC_Pericyte","41"="SMC_Pericyte",

  "39"="Cardiomyocyte",

  "8"="Neutrophil","43"="Neutrophil",

  "38"="Schwann Cell",

  "15"="B","42"="B","47"="B",

  "16"="T","23"="T","29"="T",

  "26"="NK",

  "14"="Lymph-Endo","30"="Lymph-Endo",

  "21"="Epicardium"
)

seuobj_final <- annotate_by_cluster(
  seu      = seuobj_final,
  mapping  = cluster_to_celltype,
  new_col  = "CellType1"
)

desired_levels <- c(
  "Fibroblast","Endothelium","Lymph-Endo","SMC_Pericyte","Cardiomyocyte",
  "Schwann Cell","Epicardium","Macrophage","Monocyte","cDC","pDC",
  "Neutrophil","T","NK","B","Mast Cell","Erythrocyte","Cycling Cell","Other"
)
seuobj_final$CellType1 <- factor(seuobj_final$CellType1, levels = desired_levels)

write.csv(
  table(seuobj_final$CellType1, seuobj_final$Treatment),
  file = "r_objects/celltype_by_treatment_table.csv"
)


## ---------------------------------------------------------------------------
## 11. Markers after annotation
## ---------------------------------------------------------------------------

markers_by_celltype <- find_all_markers_wrapper(
  seu   = seuobj_final,
  ident = "CellType1"
)
write.csv(
  markers_by_celltype,
  file = "r_objects/markers_by_CellType1.csv",
  row.names = FALSE
)


## ---------------------------------------------------------------------------
## 12. Visualization
## ---------------------------------------------------------------------------

p1 <- CellDimPlot(seuobj_final, group.by="CellType1", label=TRUE, reduction="UMAP") +
  plot_theme_common()
save_plot(p1, "umap_celltype1.pdf", width=10, height=8)

p2 <- CellDimPlot(seuobj_final, group.by="CellType1", label=TRUE, reduction="TSNE") +
  plot_theme_common()
save_plot(p2, "tsne_celltype1.pdf", width=10, height=8)


## ---------------------------------------------------------------------------
## 13. Save final object
## ---------------------------------------------------------------------------

saveRDS(seuobj_final, "r_objects/seuobj_final.Rds")
message("Whole-heart pipeline finished. Output saved to r_objects/seuobj_final.Rds")

###############################################################################
# End of whole_heart_main.R
###############################################################################
