###############################################################################
# cd45_neutrophil_subcluster.R
#
# Neutrophil subclustering pipeline
#
# Input:
#   r_objects/seuobj_cd45_final.Rds
#
# Steps:
#   1) Subset Neutrophils from CD45 object
#   2) Round 1: SCTransform → PCA → Harmony → clustering
#   3) Remove unwanted clusters (2,5,10,12,13,14)
#   4) Round 2: SCTransform → PCA → Harmony → clustering
#   5) Supervised annotation → NeutType1
#   6) Export markers and plots
#
# Output:
#   r_objects/seuobj_neutrophil_final.Rds
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
## 1. Load CD45 object & subset neutrophils
## ---------------------------------------------------------------------------

seu_cd45 <- readRDS("r_objects/seuobj_cd45_final.Rds")

if (!"CellType1" %in% colnames(seu_cd45@meta.data)) {
  stop("CellType1 not found. Please run cd45_main.R first.")
}

seu_neut <- subset(seu_cd45, subset = CellType1 == "Neutrophil")
seu_neut$CellType1 <- droplevels(seu_neut$CellType1)

DefaultAssay(seu_neut) <- "RNA"

message("Neutrophils detected: ", ncol(seu_neut))

## ---------------------------------------------------------------------------
## 2. Round 1: SCTransform → PCA → Harmony → clustering
## ---------------------------------------------------------------------------

dims1 <- 1:10
res1  <- 0.8

seu_neut <- run_sct_pca_harmony(
  seu             = seu_neut,
  vars_to_regress = "CC.Difference",
  npcs            = 50,
  dims_use        = dims1,
  harmony_group   = "orig.ident"
)

seu_neut <- run_clustering(
  seu        = seu_neut,
  dims_use   = dims1,
  resolution = res1
)

## ---------------------------------------------------------------------------
## 3. Global marker analysis
## ---------------------------------------------------------------------------

Idents(seu_neut) <- "seurat_clusters"

markers_global <- find_all_markers_wrapper(
  seu     = seu_neut,
  ident   = "seurat_clusters",
  min_pct = 0.25
)

write.csv(
  markers_global,
  file = "r_objects/neutrophil_markers_global.csv",
  row.names = FALSE
)

p_global <- DimPlot(
  seu_neut,
  reduction = "umap",
  label     = TRUE
) +
  ggtitle("Neutrophils: global clusters") +
  plot_theme_common()

save_plot(p_global, "neutrophil_umap_global.pdf", width = 7, height = 6)

## ---------------------------------------------------------------------------
## 4. Remove unwanted clusters
## ---------------------------------------------------------------------------

to_remove <- c(2,5,10,12,13,14)

clusters_all <- as.numeric(as.character(seu_neut$seurat_clusters))
keep_cells <- colnames(seu_neut)[!(clusters_all %in% to_remove)]

seu_neut_filt <- subset(seu_neut, cells = keep_cells)

message("Filtered neutrophils: ", ncol(seu_neut_filt))

## ---------------------------------------------------------------------------
## 5. Round 2: SCTransform → PCA → Harmony → clustering
## ---------------------------------------------------------------------------

dims2 <- 1:10
res2  <- 0.8

seu_neut_filt <- run_sct_pca_harmony(
  seu             = seu_neut_filt,
  vars_to_regress = "CC.Difference",
  npcs            = 50,
  dims_use        = dims2,
  harmony_group   = "orig.ident"
)

seu_neut_filt <- run_clustering(
  seu        = seu_neut_filt,
  dims_use   = dims2,
  resolution = res2
)

Idents(seu_neut_filt) <- "seurat_clusters"

markers_sub <- find_all_markers_wrapper(
  seu     = seu_neut_filt,
  ident   = "seurat_clusters",
  min_pct = 0.25
)

write.csv(
  markers_sub,
  file = "r_objects/neutrophil_markers_filtered_round2.csv",
  row.names = FALSE
)

p_sub <- DimPlot(
  seu_neut_filt,
  reduction = "umap",
  label     = TRUE
) +
  ggtitle("Neutrophils: filtered clusters") +
  plot_theme_common()

save_plot(p_sub, "neutrophil_umap_filtered.pdf", width = 7, height = 6)

## ---------------------------------------------------------------------------
## 6. Supervised annotation → NeutType1
## ---------------------------------------------------------------------------

neut_map <- c(
  "0"  = "Csf3r-hi Neut",
  "1"  = "Icam1+ Neut",
  "2"  = "Cd14-hi Neut",
  "3"  = "Mmp8+ Neut",
  "4"  = "Cstdc4+ Neut",
  "5"  = "Icam1+ Neut",
  "6"  = "Mmp8+ Neut",
  "7"  = "Isg15+ Neut",
  "8"  = "Icam1+ Neut",
  "9"  = "Mmp8+ Neut",
  "10" = "Mmp8+ Neut"
)

seu_neut_filt <- annotate_by_cluster(
  seu     = seu_neut_filt,
  mapping = neut_map,
  new_col = "NeutType1"
)

Idents(seu_neut_filt) <- "NeutType1"

markers_neuttype1 <- find_all_markers_wrapper(
  seu     = seu_neut_filt,
  ident   = "NeutType1",
  min_pct = 0.25
)

write.csv(
  markers_neuttype1,
  file = "r_objects/neutrophil_markers_NeutType1.csv",
  row.names = FALSE
)

p_neuttype1 <- DimPlot(
  seu_neut_filt,
  reduction = "umap",
  group.by  = "NeutType1",
  label     = TRUE
) +
  ggtitle("NeutType1 (supervised annotation)") +
  plot_theme_common()

save_plot(p_neuttype1, "neutrophil_umap_NeutType1.pdf", width = 9, height = 6)

## ---------------------------------------------------------------------------
## 7. Save final object
## ---------------------------------------------------------------------------

saveRDS(
  seu_neut_filt,
  file = "r_objects/seuobj_neutrophil_final.Rds"
)

message("Neutrophil subcluster pipeline completed. Output: r_objects/seuobj_neutrophil_final.Rds")

###############################################################################
# End of cd45_neutrophil_subcluster.R
###############################################################################
