###############################################################################
# cd45_b_subcluster.R
#
# B-cell subclustering pipeline
#
# Input :
#   r_objects/seuobj_cd45_final.Rds
#
# Steps:
#   1) Subset B cells from the CD45 object
#   2) Optional removal of cluster 30 (legacy doublet cluster)
#   3) SCTransform → PCA → Harmony → clustering
#   4) Supervised BType1 annotation
#   5) Export marker tables and plots
#
# Output:
#   r_objects/seuobj_b_final.Rds
#   r_objects/b_markers_by_cluster.csv
#   r_objects/b_markers_by_BType1.csv
###############################################################################

## ---------------------------------------------------------------------------
## 0. Packages, utils, and options
## ---------------------------------------------------------------------------

required_pkgs <- c("Seurat","harmony","dplyr","ggplot2","patchwork")
missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)

source("utils_scRNA.R")

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)

dir.create("output_plots", showWarnings = FALSE)
dir.create("r_objects", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. Load CD45 object & subset B cells
## ---------------------------------------------------------------------------

seu_cd45 <- readRDS("r_objects/seuobj_cd45_final.Rds")

if (!"CellType1" %in% colnames(seu_cd45@meta.data)) {
  stop("CellType1 not found. Please run cd45_main.R first.")
}

seu_b <- subset(seu_cd45, subset = CellType1 == "B")
seu_b$CellType1 <- droplevels(seu_b$CellType1)

DefaultAssay(seu_b) <- "RNA"

message("B cells detected: ", ncol(seu_b))

## ---------------------------------------------------------------------------
## 2. Remove cluster 30 (legacy doublet)
## ---------------------------------------------------------------------------

if ("seurat_clusters" %in% colnames(seu_b@meta.data)) {
  if ("30" %in% as.character(seu_b$seurat_clusters)) {
    seu_b <- subset(seu_b, subset = seurat_clusters != "30")
    message("Cluster 30 removed.")
  }
}

## ---------------------------------------------------------------------------
## 3. SCTransform → PCA → Harmony → clustering
## ---------------------------------------------------------------------------

dims_b <- 1:10
res_b  <- 0.8

seu_b <- run_sct_pca_harmony(
  seu             = seu_b,
  vars_to_regress = "CC.Difference",
  npcs            = 50,
  dims_use        = dims_b,
  harmony_group   = "orig.ident"
)

seu_b <- run_clustering(
  seu        = seu_b,
  dims_use   = dims_b,
  resolution = res_b
)

## ---------------------------------------------------------------------------
## 4. Marker table (cluster-level)
## ---------------------------------------------------------------------------

Idents(seu_b) <- "seurat_clusters"

markers_cluster <- find_all_markers_wrapper(
  seu     = seu_b,
  ident   = "seurat_clusters",
  min_pct = 0.25
)

write.csv(
  markers_cluster,
  file      = "r_objects/b_markers_by_cluster.csv",
  row.names = FALSE
)

p_clusters <- DimPlot(
  seu_b,
  reduction = "umap",
  label     = TRUE
) +
  ggtitle("B-cell clusters") +
  plot_theme_common()

save_plot(p_clusters, "b_umap_clusters.pdf", width = 7, height = 6)

## ---------------------------------------------------------------------------
## 5. Supervised annotation → BType1
## ---------------------------------------------------------------------------
# Mapping must be checked against your real biological definitions

b_cluster_map <- c(
  "0"  = "B-Naive",
  "1"  = "B-Naive",
  "2"  = "B-Act2",
  "3"  = "B-Act1",
  "4"  = "B-Iglc",
  "5"  = "B-Memory",
  "6"  = "B-Naive",
  "7"  = "B-Act1",
  "8"  = "B-Act2",
  "9"  = "B-Mye",
  "10" = "B-Plasma"
)

seu_b <- annotate_by_cluster(
  seu     = seu_b,
  mapping = b_cluster_map,
  new_col = "BType1"
)

Idents(seu_b) <- "BType1"

## ---------------------------------------------------------------------------
## 6. Marker table (BType1-level)
## ---------------------------------------------------------------------------

markers_btype <- find_all_markers_wrapper(
  seu     = seu_b,
  ident   = "BType1",
  min_pct = 0.25
)

write.csv(
  markers_btype,
  file      = "r_objects/b_markers_by_BType1.csv",
  row.names = FALSE
)

p_btype <- DimPlot(
  seu_b,
  reduction = "umap",
  group.by  = "BType1",
  label     = TRUE
) +
  ggtitle("B-cell subtypes (BType1)") +
  plot_theme_common()

save_plot(p_btype, "b_umap_BType1.pdf", width = 9, height = 6)

## ---------------------------------------------------------------------------
## 7. Save final object
## ---------------------------------------------------------------------------

saveRDS(seu_b, file = "r_objects/seuobj_b_final.Rds")

message("B-cell subcluster pipeline completed. Final object: r_objects/seuobj_b_final.Rds")

###############################################################################
# End of cd45_b_subcluster.R
###############################################################################
