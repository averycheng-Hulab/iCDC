###############################################################################
# whole_heart_fibroblast.R
#
# Fibroblast subclustering for whole-heart dataset
# Input:
#   r_objects/seuobj_final.Rds (from whole_heart_main.R)
#
# Steps:
#   * subset Fibroblast
#   * SCTransform (regress CC.Difference)
#   * PCA → Harmony → Clustering
#   * UMAP / TSNE visualization
#   * Manual annotation (FibType)
#   * Marker detection (cluster & FibType)
#
# Output:
#   r_objects/seuobj_fibroblast.Rds
#   fibro_markers_by_cluster.csv
#   fibro_markers_by_FibType.csv
###############################################################################

## ---------------------------------------------------------------------------
## 0. Packages & utils
## ---------------------------------------------------------------------------

required_pkgs <- c(
  "Seurat","harmony","dplyr","ggplot2","patchwork",
  "clustree","SCP"
)
missing <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse=", "))
}

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clustree)
library(SCP)

source("utils_scRNA.R")

dir.create("fibroblast_plots", showWarnings = FALSE)
dir.create("r_objects", showWarnings = FALSE)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)


## ---------------------------------------------------------------------------
## 1. Load object & subset Fibroblasts
## ---------------------------------------------------------------------------

seu_whole <- readRDS("r_objects/seuobj_final.Rds")

if (!"CellType1" %in% colnames(seu_whole@meta.data)) {
  stop("CellType1 not found. Run whole_heart_main.R first.")
}

fibro <- subset(seu_whole, subset = CellType1 == "Fibroblast")
fibro$CellType1 <- droplevels(fibro$CellType1)

message("Fibroblast cells: ", ncol(fibro))

# whole_heart_main.R already computed: S.Score / G2M.Score / Phase / CC.Difference
cc_cols <- c("S.Score","G2M.Score","Phase","CC.Difference")
miss <- setdiff(cc_cols, colnames(fibro@meta.data))
if (length(miss) > 0) stop("Missing CC metadata: ", paste(miss, collapse=", "))


## ---------------------------------------------------------------------------
## 2. SCTransform → PCA → Harmony (utils_scRNA)
## ---------------------------------------------------------------------------

pc_dims_fib <- 1:30

seu_fibro <- run_sct_pca_harmony(
  fibro,
  vars_to_regress = "CC.Difference",
  npcs            = 50,
  dims_use        = pc_dims_fib,
  harmony_group   = "orig.ident"
)


## ---------------------------------------------------------------------------
## 3. Clustering (utils_scRNA)
## ---------------------------------------------------------------------------

resolutions_fib <- c(0.4, 0.8, 1.2, 1.6)

for (res in resolutions_fib) {
  seu_fibro <- run_clustering(
    seu_fibro,
    dims_use   = pc_dims_fib,
    resolution = res
  )
}

# clustree
ct <- clustree(seu_fibro@meta.data, prefix="SCT_snn_res.")
save_plot(ct, "fibro_clustree.pdf", width=12, height=10)

# choose SCT_snn_res.0.8
if ("SCT_snn_res.0.8" %in% colnames(seu_fibro@meta.data)) {
  seu_fibro$seurat_clusters <- seu_fibro$SCT_snn_res.0.8
}


## ---------------------------------------------------------------------------
## 4. UMAP overview
## ---------------------------------------------------------------------------

save_plot(
  DimPlot(seu_fibro, reduction="umap", label=TRUE) +
    ggtitle("Fibroblast clusters") +
    plot_theme_common(),
  file.path("fibroblast_plots","fibro_umap_clusters.pdf"),
  width=8, height=6
)

save_plot(
  DimPlot(seu_fibro, reduction="umap", group.by="Treatment") +
    ggtitle("Fibroblasts by Treatment") +
    plot_theme_common(),
  file.path("fibroblast_plots","fibro_umap_treatment.pdf"),
  width=8, height=6
)


## ---------------------------------------------------------------------------
## 5. Markers (pre-annotation)
## ---------------------------------------------------------------------------

Idents(seu_fibro) <- "seurat_clusters"
markers_pre <- find_all_markers_wrapper(seu_fibro)
write.csv(markers_pre, "r_objects/fibro_markers_by_cluster.csv", row.names=FALSE)


## ---------------------------------------------------------------------------
## 6. Manual annotation (FibType)
## ---------------------------------------------------------------------------

# remove non-fibroblast contamination (endothelial-like cluster 17)
seu_fibro <- subset(seu_fibro, subset = seurat_clusters != 17)

fibro_map <- c(
  "0"="F-SL","1"="F-SL","2"="F-SL","6"="F-SL","8"="F-SL","9"="F-SL",
  "3"="F-SH","4"="F-SH","7"="F-SH","16"="F-SH",
  "10"="F-Myo","11"="F-Myo","13"="F-Myo","18"="F-Myo",
  "5"="F-Act","12"="F-Act",
  "14"="F-IFNs",
  "15"="F-IR"
)

seu_fibro <- annotate_by_cluster(
  seu_fibro,
  mapping = fibro_map,
  new_col = "FibType"
)

seu_fibro$FibType <- factor(
  seu_fibro$FibType,
  levels = c("F-SL","F-SH","F-Myo","F-Act","F-IFNs","F-IR")
)


## ---------------------------------------------------------------------------
## 7. Markers (post-annotation)
## ---------------------------------------------------------------------------

Idents(seu_fibro) <- "FibType"
markers_post <- find_all_markers_wrapper(seu_fibro, ident="FibType")
write.csv(markers_post, "r_objects/fibro_markers_by_FibType.csv", row.names=FALSE)


## ---------------------------------------------------------------------------
## 8. FibType visualization
## ---------------------------------------------------------------------------

save_plot(
  CellDimPlot(seu_fibro, group.by="FibType", reduction="UMAP",
              label=TRUE, label_insitu=TRUE, pt.size=0.1) +
    plot_theme_common(),
  file.path("fibroblast_plots","fibro_umap_FibType.pdf"),
  width=9, height=7
)

save_plot(
  CellDimPlot(seu_fibro, group.by="FibType", reduction="TSNE",
              label=TRUE, label_insitu=TRUE, pt.size=0.1) +
    plot_theme_common(),
  file.path("fibroblast_plots","fibro_tsne_FibType.pdf"),
  width=9, height=7
)

save_plot(
  CellDimPlot(seu_fibro, group.by="FibType", split.by="Treatment",
              reduction="UMAP", label=TRUE, pt.size=0.08),
  file.path("fibroblast_plots","fibro_umap_FibType_by_treatment.pdf"),
  width=14, height=8
)

save_plot(
  CellStatPlot(seu_fibro, stat.by="FibType", group.by="Treatment",
               plot_type="trend", label=TRUE) +
    plot_theme_common(),
  file.path("fibroblast_plots","fibro_cellstat_by_treatment.pdf"),
  width=10, height=6
)

save_plot(
  CellStatPlot(seu_fibro, stat.by="FibType", group.by="orig.ident",
               plot_type="trend", label=TRUE) +
    plot_theme_common(),
  file.path("fibroblast_plots","fibro_cellstat_by_sample.pdf"),
  width=14, height=8
)


## ---------------------------------------------------------------------------
## 9. Save object
## ---------------------------------------------------------------------------

saveRDS(seu_fibro, "r_objects/seuobj_fibroblast.Rds")
message("Fibroblast subclustering finished → r_objects/seuobj_fibroblast.Rds")

###############################################################################
# End of whole_heart_fibroblast.R
###############################################################################
