###############################################################################
# cd45_tnk_subcluster.R
#
# CD45+ T / NK / ILC subclustering pipeline
#
# Input :
#   r_objects/seuobj_cd45_final.Rds
#
# Steps:
#   1) Subset T and NK cells from the CD45 object
#   2) Initial T/NK clustering (dim = 20, res = 0.8) to identify/remove
#      suspected myeloid-like doublet clusters
#   3) Re-cluster cleaned T/NK cells (dim = 15, res = 1.6)
#   4) Supervised annotation → TnkType1
#   5) Fine re-clustering of T-Cyc and T-IFNa → TnkType_sel
#   6) Merge refined labels back → TnkType1_refined
#   7) Collapse NK-Cyc into NK → TnkType1_final
#   8) Save final object and marker tables
#
# Output:
#   r_objects/seuobj_tnk_final.Rds
#   r_objects/tnk_markers_by_cluster_preQC.csv
#   r_objects/tnk_markers_by_cluster_main.csv
#   r_objects/tnk_markers_by_TnkType1.csv
#   r_objects/tnk_markers_by_TnkType1_final.csv
###############################################################################

## ---------------------------------------------------------------------------
## 0. Packages, utils, options
## ---------------------------------------------------------------------------

required_pkgs <- c(
  "Seurat",
  "harmony",
  "dplyr",
  "ggplot2",
  "patchwork"
)

missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
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

dir.create("r_objects", showWarnings = FALSE)
dir.create("output_plots", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. Load CD45 object and subset T / NK
## ---------------------------------------------------------------------------

seuobj_cd45_final <- readRDS("r_objects/seuobj_cd45_final.Rds")

if (!"CellType1" %in% colnames(seuobj_cd45_final@meta.data)) {
  stop("CellType1 not found in metadata. Please run cd45_main.R first.")
}

seu_tnk <- subset(
  seuobj_cd45_final,
  subset = CellType1 %in% c("T", "NK")
)
seu_tnk$CellType1 <- droplevels(seu_tnk$CellType1)

message("T/NK cells (raw subset): ", ncol(seu_tnk))

DefaultAssay(seu_tnk) <- "RNA"

## ---------------------------------------------------------------------------
## 2. Initial T/NK clustering (for doublet removal)
##    dim = 20, res = 0.8
## ---------------------------------------------------------------------------

pc_dims_pre <- 1:20
res_pre     <- 0.8

seu_tnk_pre <- run_sct_pca_harmony(
  seu             = seu_tnk,
  vars_to_regress = "CC.Difference",
  npcs            = 50,
  dims_use        = pc_dims_pre,
  harmony_group   = "orig.ident"
)

seu_tnk_pre <- run_clustering(
  seu        = seu_tnk_pre,
  dims_use   = pc_dims_pre,
  resolution = res_pre
)

markers_pre <- find_all_markers_wrapper(
  seu     = seu_tnk_pre,
  ident   = "seurat_clusters",
  min_pct = 0.25
)
write.csv(
  markers_pre,
  file      = "r_objects/tnk_markers_by_cluster_preQC.csv",
  row.names = FALSE
)

p_pre <- DimPlot(
  seu_tnk_pre,
  reduction = "umap",
  label     = TRUE
) +
  ggtitle("T/NK pre-clusters (doublet screening)") +
  plot_theme_common()
save_plot(p_pre, "tnk_umap_preclusters.pdf", width = 8, height = 6)

## ---------------------------------------------------------------------------
## 3. Remove suspected myeloid-like doublet clusters
## ---------------------------------------------------------------------------

doublet_clusters <- c("7", "10", "14", "17")

seu_tnk_clean <- subset(
  seu_tnk_pre,
  subset = !(seurat_clusters %in% doublet_clusters)
)

message("T/NK cells after removing suspected doublet clusters: ", ncol(seu_tnk_clean))

## ---------------------------------------------------------------------------
## 4. Main T/NK subclustering
##    dim = 15, res = 1.6
## ---------------------------------------------------------------------------

pc_dims_main <- 1:15
res_main     <- 1.6

DefaultAssay(seu_tnk_clean) <- "RNA"

seu_tnk_main <- run_sct_pca_harmony(
  seu             = seu_tnk_clean,
  vars_to_regress = "CC.Difference",
  npcs            = 50,
  dims_use        = pc_dims_main,
  harmony_group   = "orig.ident"
)

seu_tnk_main <- run_clustering(
  seu        = seu_tnk_main,
  dims_use   = pc_dims_main,
  resolution = res_main
)

markers_main <- find_all_markers_wrapper(
  seu     = seu_tnk_main,
  ident   = "seurat_clusters",
  min_pct = 0.25
)
write.csv(
  markers_main,
  file      = "r_objects/tnk_markers_by_cluster_main.csv",
  row.names = FALSE
)

p_main <- DimPlot(
  seu_tnk_main,
  reduction = "umap",
  label     = TRUE
) +
  ggtitle("T/NK main clusters (dim = 15, res = 1.6)") +
  plot_theme_common()
save_plot(p_main, "tnk_umap_clusters_main.pdf", width = 8, height = 6)

## ---------------------------------------------------------------------------
## 5. Supervised annotation → TnkType1
## ---------------------------------------------------------------------------

cluster_to_TnkType1 <- c(
  "0"  = "T-Naive",
  "1"  = "T-Naive",
  "2"  = "Tgd",
  "3"  = "T-Helper",
  "4"  = "NK",
  "5"  = "NKT",
  "6"  = "NK",
  "7"  = "ILC",
  "8"  = "T-Toxic",
  "9"  = "T-Toxic",
  "10" = "T-DN",
  "11" = "T-Cyc",
  "12" = "T-Reg",
  "13" = "T-Toxic",
  "14" = "T-IFNa",
  "15" = "Tgd",
  "16" = "NK",
  "17" = "NK-Cyc",
  "18" = "T-CM",
  "19" = "Tgd",
  "20" = "Other",
  "21" = "T-Toxic",
  "22" = "Tgd",
  "23" = "T-CM"
)

seu_tnk_main <- annotate_by_cluster(
  seu     = seu_tnk_main,
  mapping = cluster_to_TnkType1,
  new_col = "TnkType1"
)

tnk_levels <- c(
  "T-Helper","T-Toxic","Tgd","T-Reg","T-Naive","T-CM",
  "T-DN","T-IFNa","T-Cyc","NKT","NK","NK-Cyc","ILC","Other"
)
seu_tnk_main$TnkType1 <- factor(seu_tnk_main$TnkType1, levels = tnk_levels)

markers_TnkType1 <- find_all_markers_wrapper(
  seu     = seu_tnk_main,
  ident   = "TnkType1",
  min_pct = 0.25
)
write.csv(
  markers_TnkType1,
  file      = "r_objects/tnk_markers_by_TnkType1.csv",
  row.names = FALSE
)

p_tnk_type1 <- DimPlot(
  seu_tnk_main,
  reduction = "umap",
  group.by  = "TnkType1",
  label     = TRUE
) +
  ggtitle("T/NK subtypes (TnkType1, pre-refine)") +
  plot_theme_common()
save_plot(p_tnk_type1, "tnk_umap_TnkType1_pre_refine.pdf", width = 9, height = 7)

## ---------------------------------------------------------------------------
## 6. Fine re-clustering of T-Cyc and T-IFNa → TnkType_sel
## ---------------------------------------------------------------------------

seu_tsel <- subset(
  seu_tnk_main,
  subset = TnkType1 %in% c("T-Cyc", "T-IFNa")
)

message("Cells in T-Cyc / T-IFNa for fine re-clustering: ", ncol(seu_tsel))

if (ncol(seu_tsel) > 0) {
  DefaultAssay(seu_tsel) <- "RNA"

  pc_dims_sel <- 1:10
  res_sel     <- 0.8

  seu_tsel <- run_sct_pca_harmony(
    seu             = seu_tsel,
    vars_to_regress = "CC.Difference",
    npcs            = 50,
    dims_use        = pc_dims_sel,
    harmony_group   = "orig.ident"
  )

  seu_tsel <- run_clustering(
    seu        = seu_tsel,
    dims_use   = pc_dims_sel,
    resolution = res_sel
  )

  cluster_to_TnkType_sel <- c(
    "1" = "T-Helper",
    "2" = "T-Helper",
    "0" = "T-Toxic",
    "3" = "T-Toxic",
    "4" = "T-Toxic",
    "5" = "T-Toxic",
    "8" = "Tgd",
    "6" = "NKT",
    "7" = "T-DN"
  )

  seu_tsel <- annotate_by_cluster(
    seu     = seu_tsel,
    mapping = cluster_to_TnkType_sel,
    new_col = "TnkType_sel"
  )

  sel_levels <- c("T-Helper","T-Toxic","Tgd","NKT","T-DN","Other")
  seu_tsel$TnkType_sel <- factor(seu_tsel$TnkType_sel, levels = sel_levels)

  seu_tnk_main <- merge_annotation_back(
    parent_obj  = seu_tnk_main,
    refined_obj = seu_tsel,
    parent_col  = "TnkType1",
    refined_col = "TnkType_sel",
    new_col     = "TnkType1_refined"
  )
} else {
  seu_tnk_main$TnkType1_refined <- seu_tnk_main$TnkType1
}

## ---------------------------------------------------------------------------
## 7. Collapse NK-Cyc into NK → TnkType1_final
## ---------------------------------------------------------------------------

base_type <- as.character(seu_tnk_main$TnkType1_refined)
base_type[base_type == "NK-Cyc"] <- "NK"

final_levels <- c(
  "T-Helper","T-Toxic","Tgd","T-Reg","T-Naive","T-CM",
  "T-DN","NKT","NK","ILC","Other"
)

seu_tnk_main$TnkType1_final <- factor(base_type, levels = final_levels)
Idents(seu_tnk_main) <- "TnkType1_final"

markers_TnkType1_final <- find_all_markers_wrapper(
  seu     = seu_tnk_main,
  ident   = "TnkType1_final",
  min_pct = 0.25
)
write.csv(
  markers_TnkType1_final,
  file      = "r_objects/tnk_markers_by_TnkType1_final.csv",
  row.names = FALSE
)

p_tnk_final <- DimPlot(
  seu_tnk_main,
  reduction = "umap",
  group.by  = "TnkType1_final",
  label     = TRUE
) +
  ggtitle("T/NK subtypes (TnkType1_final)") +
  plot_theme_common()
save_plot(p_tnk_final, "tnk_umap_TnkType1_final.pdf", width = 9, height = 7)

## ---------------------------------------------------------------------------
## 8. Save final T/NK object
## ---------------------------------------------------------------------------

seuobj_tnk_final <- seu_tnk_main
saveRDS(seuobj_tnk_final, file = "r_objects/seuobj_tnk_final.Rds")

message("CD45+ T/NK subcluster pipeline finished. Final object: r_objects/seuobj_tnk_final.Rds")

###############################################################################
# End of cd45_tnk_subcluster.R
###############################################################################
