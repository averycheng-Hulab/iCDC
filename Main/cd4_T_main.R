###############################################################################
# cd4_T_main.R
#
# CD4+ T cell scRNA-seq pipeline (10X, Seurat + DoubletFinder + Harmony)
#
# Steps:
#   - Per-sample QC → filtering → DoubletFinder
#   - Merge → remove mito/ribo genes
#   - Cell cycle scoring
#   - SCTransform → PCA → Harmony integration
#   - Clustering + UMAP/TSNE
#   - Manual annotation (CD4 T-cell subsets)
#
# Outputs:
#   r_objects/seuobj_cd4_final.Rds
#   cd4_markers_by_cluster.csv
#   cd4_markers_by_CellType1.csv
###############################################################################

## ---------------------------------------------------------------------------
## 0. Load packages and utilities
## ---------------------------------------------------------------------------

required_pkgs <- c(
  "Seurat","DoubletFinder","harmony",
  "ggplot2","patchwork","dplyr","stringr","clustree"
)

missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse=", "))
}

library(Seurat)
library(DoubletFinder)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(clustree)

source("utils_scRNA.R")   # ★ 使用统一的 QC / CC scoring / SCT / Harmony / clustering / plots

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)

dir.create("QC", showWarnings = FALSE)
dir.create("output_plots", showWarnings = FALSE)
dir.create("r_objects", showWarnings = FALSE)


## ---------------------------------------------------------------------------
## 1. Sample names & paths
## ---------------------------------------------------------------------------

base_dir <- "TCR_mRNA"
samples  <- list.files(base_dir)

treatment_map <- setNames(c("iCDC","iCDC","MI","MI"), samples)
organ_map     <- setNames(c("Blood","Heart","Blood","Heart"), samples)

matrix_dirs <- file.path(base_dir, samples, "sample_filtered_feature_bc_matrix")
stopifnot(length(samples) == length(matrix_dirs))


## ---------------------------------------------------------------------------
## 2. Per-sample QC + DoubletFinder  → seuobj_list
## ---------------------------------------------------------------------------

seuobj_list <- vector("list", length(samples))
names(seuobj_list) <- samples

for (i in seq_along(samples)) {
  seuobj_list[[i]] <- process_sample_scRNA(
    matrix_dir      = matrix_dirs[i],
    sample_name     = samples[i],
    treatment_label = treatment_map[[samples[i]]],
    min_features    = 300,     # utils default
    min_counts      = 500,     # utils default
    max_mito        = 5,       # utils default
    hb_thr          = 5,       # utils default
    pc_dims         = 1:40,
    doublet_rate    = 0.008
  )
  seuobj_list[[i]]$Organ <- organ_map[[samples[i]]]
}


## ---------------------------------------------------------------------------
## 3. Merge → seuobj_raw
## ---------------------------------------------------------------------------

seuobj_merged <- seuobj_list[[1]]
for (i in 2:length(seuobj_list)) {
  seuobj_merged <- merge(seuobj_merged, y = seuobj_list[[i]])
}

BasicData <- GetAssayData(seuobj_merged, assay="RNA", slot="counts")
meta_min  <- seuobj_merged@meta.data %>% dplyr::select(orig.ident,Treatment,Organ)

seuobj_raw <- CreateSeuratObject(counts=BasicData, meta.data=meta_min)

seuobj_raw$orig.ident <- factor(seuobj_raw$orig.ident, levels=samples)
seuobj_raw$Treatment  <- factor(seuobj_raw$Treatment, levels=c("iCDC","MI"))
seuobj_raw$Organ      <- factor(seuobj_raw$Organ, levels=c("Blood","Heart"))


## ---------------------------------------------------------------------------
## 4. Remove mito + ribosomal genes
## ---------------------------------------------------------------------------

seuobj_qc <- filter_mito_ribo(seuobj_raw)


## ---------------------------------------------------------------------------
## 5. Cell cycle scoring（utils_scRNA 封装）
## ---------------------------------------------------------------------------

data("cc.genes.updated.2019", package="Seurat")
s_genes   <- str_to_title(cc.genes.updated.2019$s.genes)
g2m_genes <- str_to_title(cc.genes.updated.2019$g2m.genes)

seuobj_qc <- compute_cell_cycle_scores(
  seuobj_qc,
  s_genes   = s_genes,
  g2m_genes = g2m_genes
)


## ---------------------------------------------------------------------------
## 6. SCTransform → PCA → Harmony（utils_scRNA 封装）
## ---------------------------------------------------------------------------

seuobj_final <- run_sct_pca_harmony(
  seuobj_qc,
  vars_to_regress = "CC.Difference",
  npcs            = 50,
  dims_use        = 1:40,
  harmony_group   = "orig.ident"
)


## ---------------------------------------------------------------------------
## 7. Clustering + TSNE + UMAP（utils_scRNA 封装）
## ---------------------------------------------------------------------------

cluster_resolutions <- seq(0.2, 1.6, 0.2)

for (res in cluster_resolutions) {
  seuobj_final <- run_clustering(
    seuobj_final,
    dims_use   = 1:40,
    resolution = res
  )
}

# clustree overview
ctree <- clustree(seuobj_final@meta.data, prefix="SCT_snn_res.")
save_plot(ctree, "cd4_clustree.pdf", width=12, height=8)

# choose resolution 0.8
seuobj_final$seurat_clusters <- seuobj_final$SCT_snn_res.0.8


## ---------------------------------------------------------------------------
## 8. Pre-annotation markers
## ---------------------------------------------------------------------------

markers_by_cluster <- find_all_markers_wrapper(seuobj_final, ident="seurat_clusters")
write.csv(markers_by_cluster, "r_objects/cd4_markers_by_cluster.csv", row.names=FALSE)

save_plot(
  DimPlot(seuobj_final, reduction="umap", label=TRUE) + plot_theme_common(),
  "cd4_umap_clusters.pdf", 7,6
)
save_plot(
  DimPlot(seuobj_final, reduction="umap", group.by="orig.ident") + plot_theme_common(),
  "cd4_umap_by_sample.pdf", 7,6
)
save_plot(
  DimPlot(seuobj_final, reduction="umap", group.by="Treatment") + plot_theme_common(),
  "cd4_umap_by_treatment.pdf", 7,6
)
save_plot(
  DimPlot(seuobj_final, reduction="umap", group.by="Organ") + plot_theme_common(),
  "cd4_umap_by_organ.pdf", 7,6
)


## ---------------------------------------------------------------------------
## 9. Manual annotation（CD4 T subsets）
## ---------------------------------------------------------------------------

cluster_to_celltype <- c(
  "0"  = "T-CM/Naive",
  "1"  = "T-Reg",
  "2"  = "T-CM/Naive",
  "3"  = "T-CM/Naive",
  "4"  = "T-CM/Naive",
  "5"  = "T-EarlyAct",
  "6"  = "T-Isg",
  "7"  = "T-CM/Naive",
  "8"  = "T-Reg",
  "9"  = "T-CM/Naive",
  "10" = "T-CM/Naive",
  "11" = "T-Trbv1",
  "12" = "T-EM",
  "13" = "T-Cycling",
  "14" = "T-Reg"
)

seuobj_final <- annotate_by_cluster(
  seuobj_final,
  mapping = cluster_to_celltype,
  new_col = "CellType1"
)

desired_levels <- c(
  "T-CM/Naive","T-Reg","T-EarlyAct","T-Isg",
  "T-Trbv1","T-EM","T-Cycling","Other"
)
seuobj_final$CellType1 <- factor(seuobj_final$CellType1, levels=desired_levels)


## ---------------------------------------------------------------------------
## 10. Post-annotation markers
## ---------------------------------------------------------------------------

markers_by_celltype <- find_all_markers_wrapper(seuobj_final, ident="CellType1")
write.csv(markers_by_celltype, "r_objects/cd4_markers_by_CellType1.csv", row.names=FALSE)


## ---------------------------------------------------------------------------
## 11. Plots by CellType1
## ---------------------------------------------------------------------------

save_plot(
  DimPlot(seuobj_final, group.by="CellType1", label=TRUE, reduction="UMAP") +
    plot_theme_common(),
  "cd4_umap_CellType1.pdf", 8,6
)
save_plot(
  DimPlot(seuobj_final, group.by="CellType1", label=TRUE, reduction="TSNE") +
    plot_theme_common(),
  "cd4_tsne_CellType1.pdf", 8,6
)


## ---------------------------------------------------------------------------
## 12. Save final Seurat object
## ---------------------------------------------------------------------------

saveRDS(seuobj_final, "r_objects/seuobj_cd4_final.Rds")
message("CD4 T pipeline finished. Output → r_objects/seuobj_cd4_final.Rds")

###############################################################################
# End of cd4_T_main.R
###############################################################################
