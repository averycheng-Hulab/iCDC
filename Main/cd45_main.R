###############################################################################
# CD45+ scRNA-seq pipeline (10X, Seurat + DoubletFinder + Harmony)
#
# - Input: 11 samples under Summary/<sample>/
# - Per sample:
#       QC → filtering → DoubletFinder → singlet-only Seurat object
# - After merge:
#       remove mito/ribo → cell cycle scoring
#       SCTransform(regress CC.Difference) → PCA → Harmony
#       clustering → UMAP / TSNE
#       manual annotation
#
# - Outputs:
#       r_objects/seuobj_final.Rds
#       markers_by_cluster.csv
#       markers_by_CellType1.csv
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
}

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
## 1. Sample names & 10X paths
## ---------------------------------------------------------------------------

samples <- c(
  "DC_1","DC_2","DC_3","DC_4",
  "IR_1","IR_2","IR_3","IR_4",
  "Sham_1","Sham_2","Sham_3"
)

treatment_map <- c(
  "DC_1"="DC","DC_2"="DC","DC_3"="DC","DC_4"="DC",
  "IR_1"="IR","IR_2"="IR","IR_3"="IR","IR_4"="IR",
  "Sham_1"="Sham","Sham_2"="Sham","Sham_3"="Sham"
)

base_dir <- "Summary"
matrix_dirs <- file.path(base_dir, samples)
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
## 4. Merge → seuobj_raw
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
seuobj_raw$Treatment  <- factor(seuobj_raw$Treatment, levels = c("Sham","IR","DC"))

## ---------------------------------------------------------------------------
## 5. Remove mito/ribo genes
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
## 7. SCTransform + PCA + Harmony
## ---------------------------------------------------------------------------

seuobj_final <- run_sct_pca_harmony(
  seu             = seuobj_qc,
  vars_to_regress = "CC.Difference",
  npcs            = 50,
  dims_use        = pc_dims,
  harmony_group   = "orig.ident"
)

## ---------------------------------------------------------------------------
## 8. Clustering + UMAP + TSNE
## ---------------------------------------------------------------------------

for (res in cluster_res) {
  seuobj_final <- run_clustering(
    seu        = seuobj_final,
    dims_use   = pc_dims,
    resolution = res
  )
}

ctree <- clustree(seuobj_final@meta.data, prefix="SCT_snn_res.")
save_plot(ctree, "cd45_clustree.pdf", width=15, height=10)

save_plot(
  DimPlot(seuobj_final, reduction="umap", label=TRUE) +
    ggtitle("CD45+: clusters") + plot_theme_common(),
  "umap_clusters.pdf"
)

save_plot(
  DimPlot(seuobj_final, reduction="umap", group.by="orig.ident") +
    ggtitle("CD45+: by sample") + plot_theme_common(),
  "umap_by_sample.pdf"
)

save_plot(
  DimPlot(seuobj_final, reduction="umap", group.by="Treatment") +
    ggtitle("CD45+: by treatment") + plot_theme_common(),
  "umap_by_treatment.pdf"
)

## ---------------------------------------------------------------------------
## 9. Pre-annotation markers (seurat_clusters)
## ---------------------------------------------------------------------------

markers_by_cluster <- find_all_markers_wrapper(
  seu   = seuobj_final,
  ident = "seurat_clusters"
)
write.csv(markers_by_cluster, "r_objects/markers_by_cluster.csv", row.names=FALSE)

## ---------------------------------------------------------------------------
## 10. Manual annotation: CellType1
## ---------------------------------------------------------------------------

cluster_to_celltype <- c(
  "0"="Mac_Mono","1"="Mac_Mono","2"="B","3"="Mac_Mono","4"="Mac_Mono",
  "5"="Mac_Mono","6"="B","7"="Mac_Mono","8"="Mac_Mono","9"="Mac_Mono",
  "10"="Mac_Mono","11"="Mac_Mono","12"="T","13"="T","14"="Neutrophil",
  "15"="Neutrophil","16"="Mac_Mono","17"="DC","18"="T","19"="Neutrophil",
  "20"="Mac_Mono","21"="Mac_Mono","22"="Neutrophil","23"="NK",
  "24"="T","25"="Mac_Mono","26"="NK","27"="DC","28"="Endothelium",
  "29"="Mac_Mono","30"="B","31"="Neutrophil","32"="T","33"="Mac_Mono",
  "34"="DC","35"="T","36"="T","37"="Endothelium","38"="DC",
  "39"="Fibroblast","40"="Other","41"="Other","42"="Other",
  "43"="Endothelium","44"="Mast","45"="Neutrophil","46"="B","47"="Mac_Mono"
)

seuobj_final <- annotate_by_cluster(
  seu     = seuobj_final,
  mapping = cluster_to_celltype,
  new_col = "CellType1"
)

desired_levels <- c(
  "Mac_Mono","DC","Neutrophil","B","T","NK",
  "Endothelium","Fibroblast","Mast","Other"
)
seuobj_final$CellType1 <- factor(seuobj_final$CellType1, levels = desired_levels)

## ---------------------------------------------------------------------------
## 11. Post-annotation markers (CellType1)
## ---------------------------------------------------------------------------

markers_by_celltype <- find_all_markers_wrapper(
  seu   = seuobj_final,
  ident = "CellType1"
)
write.csv(markers_by_celltype, "r_objects/markers_by_CellType1.csv", row.names=FALSE)

## ---------------------------------------------------------------------------
## 12. Visualization
## ---------------------------------------------------------------------------

save_plot(
  CellDimPlot(seuobj_final, group.by="CellType1", label=TRUE, reduction="UMAP") +
    plot_theme_common(),
  "umap_celltype1.pdf", width=10, height=8
)

save_plot(
  CellDimPlot(seuobj_final, group.by="CellType1", label=TRUE, reduction="TSNE") +
    plot_theme_common(),
  "tsne_celltype1.pdf", width=10, height=8
)

save_plot(
  CellStatPlot(seuobj_final, stat.by="CellType1", group.by="Treatment", plot_type="trend") +
    plot_theme_common(),
  "cellstat_by_treatment.pdf", width=10, height=6
)

save_plot(
  CellStatPlot(seuobj_final, stat.by="CellType1", group.by="orig.ident", plot_type="trend") +
    plot_theme_common(),
  "cellstat_by_sample.pdf", width=14, height=8
)

## ---------------------------------------------------------------------------
## 13. Save final object
## ---------------------------------------------------------------------------

saveRDS(seuobj_final, "r_objects/seuobj_final.Rds")
message("CD45+ pipeline finished. Output saved to r_objects/seuobj_cd45_final.Rds")

###############################################################################
# End of cd45_main.R
###############################################################################
