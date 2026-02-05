###############################################################################
# cd45_MacMonoDc_subcluster.R
#
# - Input : r_objects/seuobj_cd45_final.Rds
# - Output: r_objects/seuobj_macmonodc_final.Rds
#
# - All SCP/SingleR removed
# - All Cell composition plots removed
# - All supervised mappings preserved
# - Two mapping sources preserved (MacType2 + fine re-cluster MacType_sel)
###############################################################################

## ---------------------------------------------------------------------------
## 0. Packages, utils, options
## ---------------------------------------------------------------------------

required_pkgs <- c(
  "Seurat","harmony","dplyr","ggplot2","patchwork","stringr","RColorBrewer"
)

missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse=", "))
}

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(RColorBrewer)

source("utils_scRNA.R")

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)

dir.create("output_plots", showWarnings = FALSE)
dir.create("r_objects", showWarnings = FALSE)


## ---------------------------------------------------------------------------
## 1. Load CD45 object & subset Mac_Mono / DC
## ---------------------------------------------------------------------------

seuobj_cd45_final <- readRDS("r_objects/seuobj_cd45_final.Rds")

seuobj_macmonodc <- subset(
  seuobj_cd45_final,
  subset = CellType1 %in% c("Mac_Mono","DC")
)

seuobj_macmonodc$CellType1 <- factor(
  seuobj_macmonodc$CellType1,
  levels = c("Mac_Mono","DC")
)

message("Mac/Mono/DC cells: ", ncol(seuobj_macmonodc))


## ---------------------------------------------------------------------------
## 2. Global SCTransform + PCA + Harmony + clustering
## ---------------------------------------------------------------------------

pc_dims_global <- 1:40
res_global     <- 2.0

# CC.Difference 已在 cd45_main.R 中计算好，直接回归
seuobj_macmonodc <- run_sct_pca_harmony(
  seu   = seuobj_macmonodc,
  vars_to_regress = "CC.Difference",
  npcs  = 50,
  dims_use = pc_dims_global,
  harmony_group = "orig.ident"
)

seuobj_macmonodc <- run_clustering(
  seu       = seuobj_macmonodc,
  dims_use  = pc_dims_global,
  resolution = res_global
)

## ---------------------------------------------------------------------------
## 3. Markers by cluster（全局 Mac/Mono/DC）
## ---------------------------------------------------------------------------

markers_cluster <- find_all_markers_wrapper(
  seu   = seuobj_macmonodc,
  ident = "seurat_clusters",
  min_pct = 0.25
)
write.csv(
  markers_cluster,
  file      = file.path("r_objects", "macmonodc_markers_by_cluster.csv"),
  row.names = FALSE
)

save_plot(
  DimPlot(seuobj_macmonodc, reduction="umap", label=TRUE) +
    ggtitle("Mac/Mono/DC Clusters") + plot_theme_common(),
  filename = file.path("output_plots", "macmonodc_umap_clusters.pdf"),
  width = 8,
  height = 6
)


## ---------------------------------------------------------------------------
## 4. Supervised MacType1 mapping
## ---------------------------------------------------------------------------

mac_cluster_map <- c(
  "14"="Lyve1+ Mac","26"="Lyve1+ Mac","15"="Lyve1+ Mac","16"="Lyve1+ Mac",
  "0"="Lyve1+ Mac","23"="Lyve1+ Mac","30"="Lyve1+ Mac","6"="Lyve1+ Mac",
  "2"="Lyve1+ Mac","5"="Lyve1+ Mac",
  "10"="MHCII-hi Mac","22"="MHCII-hi Mac","25"="MHCII-hi Mac",
  "4"="MHCII-hi Mac","8"="MHCII-hi Mac","12"="MHCII-hi Mac","7"="MHCII-hi Mac",
  "11"="MHCII-hi Mac","3"="MHCII-hi Mac","27"="MHCII-hi Mac","31"="MHCII-hi Mac",
  "32"="MHCII-hi Mac","1"="MHCII-hi Mac","33"="MHCII-hi Mac","29"="MHCII-hi Mac",
  "13"="Ccl8+ Mac",
  "18"="Spp1+ Mac","20"="Spp1+ Mac",
  "21"="IFN Mac",
  "28"="Fn1+ Mac",
  "36"="Saa3+ Mac",
  "19"="Prolif Mac",
  "17"="Ly6c-hi Mono",
  "37"="Ly6c-lo Mono",
  "24"="cDC1","41"="cDC1",
  "9"="cDC2",
  "34"="pDC","42"="pDC",
  "35"="Migratory DC",
  "38"="Other","39"="Other","40"="Other","43"="Other"
)

seuobj_macmonodc <- annotate_by_cluster(
  seu     = seuobj_macmonodc,
  mapping = mac_cluster_map,
  new_col = "MacType1"
)

markers_MacType1 <- find_all_markers_wrapper(
  seu   = seuobj_macmonodc,
  ident = "MacType1",
  min_pct = 0.25
)
write.csv(
  markers_MacType1,
  file      = file.path("r_objects", "macmonodc_markers_by_MacType1.csv"),
  row.names = FALSE
)

save_plot(
  DimPlot(seuobj_macmonodc, reduction="umap",
          group.by="MacType1", label=TRUE) +
    ggtitle("MacType1") + plot_theme_common(),
  filename = file.path("output_plots", "macmonodc_umap_MacType1.pdf"),
  width  = 10,
  height = 7
)


## ---------------------------------------------------------------------------
## 5. MacType2 (Mac / Mono / DC / Other)
## ---------------------------------------------------------------------------

cluster_ids <- as.numeric(as.character(seuobj_macmonodc$seurat_clusters))

index_mono  <- c(17, 37)
index_dc    <- c(24, 41, 9, 35, 34, 42)
index_other <- c(38:40, 43)
index_mac   <- setdiff(0:43, c(index_mono, index_dc, index_other))

MacType2_vec <- ifelse(
  cluster_ids %in% index_mac,   "Mac",
  ifelse(
    cluster_ids %in% index_mono, "Mono",
    ifelse(cluster_ids %in% index_dc, "DC", "Other")
  )
)

seuobj_macmonodc$MacType2 <- factor(MacType2_vec,
                                    levels = c("Mac","Mono","DC","Other"))

markers_MacType2 <- find_all_markers_wrapper(
  seu   = seuobj_macmonodc,
  ident = "MacType2",
  min_pct = 0.25
)
write.csv(
  markers_MacType2,
  file      = file.path("r_objects", "macmonodc_markers_by_MacType2.csv"),
  row.names = FALSE
)

save_plot(
  DimPlot(seuobj_macmonodc, reduction="umap",
          group.by="MacType2", label=TRUE) +
    ggtitle("MacType2 (Mac/Mono/DC/Other)") + plot_theme_common(),
  filename = file.path("output_plots", "macmonodc_umap_MacType2.pdf"),
  width  = 9,
  height = 6
)


## ---------------------------------------------------------------------------
## 6. Fine re-cluster selected Mac clusters
## ---------------------------------------------------------------------------

selected_clusters <- c("14","26","15","29","0","5","22","25","10","2")

seuobj_macmonodc_sel <- subset(
  seuobj_macmonodc,
  subset = seurat_clusters %in% selected_clusters
)

pc_dims_sel <- 1:10
res_sel     <- 1.2

seuobj_macmonodc_sel <- run_sct_pca_harmony(
  seu            = seuobj_macmonodc_sel,
  vars_to_regress = "CC.Difference",
  npcs           = 50,
  dims_use       = pc_dims_sel,
  harmony_group  = "orig.ident"
)

seuobj_macmonodc_sel <- run_clustering(
  seu       = seuobj_macmonodc_sel,
  dims_use  = pc_dims_sel,
  resolution = res_sel
)

markers_sel_cluster <- find_all_markers_wrapper(
  seu   = seuobj_macmonodc_sel,
  ident = "seurat_clusters",
  min_pct = 0.25
)
write.csv(
  markers_sel_cluster,
  file      = file.path("r_objects", "macmonodc_sel_markers_by_cluster.csv"),
  row.names = FALSE
)

save_plot(
  DimPlot(seuobj_macmonodc_sel, reduction="umap",
          group.by="seurat_clusters", label=TRUE) +
    ggtitle("Mac Selected Clusters (Re-cluster)") + plot_theme_common(),
  filename = file.path("output_plots", "macmonodc_sel_umap_clusters.pdf"),
  width  = 7,
  height = 7
)


## ---------------------------------------------------------------------------
## 7. Supervised annotation MacType_sel
## ---------------------------------------------------------------------------

mac_sel_map <- c(
  "0"="MHCII-hi Mac","5"="MHCII-hi Mac","9"="MHCII-hi Mac","11"="MHCII-hi Mac",
  "14"="MHCII-hi Mac","15"="MHCII-hi Mac","16"="MHCII-hi Mac","17"="MHCII-hi Mac",
  "18"="MHCII-hi Mac","19"="MHCII-hi Mac","1"="MHCII-hi Mac","10"="MHCII-hi Mac",
  "2"="Lyve1+ Mac","3"="Lyve1+ Mac","4"="Lyve1+ Mac","6"="Lyve1+ Mac",
  "7"="Lyve1+ Mac","8"="Lyve1+ Mac","12"="Lyve1+ Mac","13"="Lyve1+ Mac","20"="Lyve1+ Mac"
)

seuobj_macmonodc_sel <- annotate_by_cluster(
  seu     = seuobj_macmonodc_sel,
  mapping = mac_sel_map,
  new_col = "MacType_sel"
)

markers_MacType_sel <- find_all_markers_wrapper(
  seu   = seuobj_macmonodc_sel,
  ident = "MacType_sel",
  min_pct = 0.25
)
write.csv(
  markers_MacType_sel,
  file      = file.path("r_objects", "macmonodc_sel_markers_by_MacType_sel.csv"),
  row.names = FALSE
)

save_plot(
  DimPlot(seuobj_macmonodc_sel, reduction="umap",
          group.by="MacType_sel", label=TRUE) +
    ggtitle("MacType_sel") + plot_theme_common(),
  filename = file.path("output_plots", "macmonodc_sel_umap_MacType_sel.pdf"),
  width  = 7,
  height = 7
)


## ---------------------------------------------------------------------------
## 8. Merge MacType_sel back → MacType1_refined
## ---------------------------------------------------------------------------

seuobj_macmonodc <- merge_annotation_back(
  parent_obj = seuobj_macmonodc,
  refined_obj = seuobj_macmonodc_sel,
  parent_col = "MacType1",
  refined_col = "MacType_sel",
  new_col    = "MacType1_refined"
)

# 设定 MacType1_refined 的 factor 顺序
seuobj_macmonodc$MacType1_refined <- factor(
  seuobj_macmonodc$MacType1_refined,
  levels = c(
    "Lyve1+ Mac","MHCII-hi Mac","Ccl8+ Mac","Spp1+ Mac",
    "IFN Mac","Fn1+ Mac","Saa3+ Mac","Prolif Mac",
    "Ly6c-hi Mono","Ly6c-lo Mono",
    "cDC1","cDC2","pDC","Migratory DC","Other"
  )
)

markers_MacType1_refined <- find_all_markers_wrapper(
  seu   = seuobj_macmonodc,
  ident = "MacType1_refined",
  min_pct = 0.25
)
write.csv(
  markers_MacType1_refined,
  file      = file.path("r_objects", "macmonodc_markers_by_MacType1_refined.csv"),
  row.names = FALSE
)

save_plot(
  DimPlot(seuobj_macmonodc, reduction="umap",
          group.by="MacType1_refined", label=TRUE) +
    ggtitle("MacType1_refined") + plot_theme_common(),
  filename = file.path("output_plots", "macmonodc_umap_MacType1_refined.pdf"),
  width  = 10,
  height = 7
)


## ---------------------------------------------------------------------------
## 9. Add Ccr2-based split → MacType1_final
## ---------------------------------------------------------------------------

if ("Ccr2" %in% rownames(seuobj_macmonodc)) {
  ccr2_expr <- FetchData(seuobj_macmonodc, vars="Ccr2")[,1]
  seuobj_macmonodc$Ccr2_exp <- ccr2_expr
  seuobj_macmonodc$Ccr2_pos <- factor(
    ifelse(ccr2_expr > 0, "TRUE", "FALSE"),
    levels = c("TRUE","FALSE")
  )

  base_type <- as.character(seuobj_macmonodc$MacType1_refined)
  is_mhcii  <- base_type == "MHCII-hi Mac"
  is_ccr2p  <- seuobj_macmonodc$Ccr2_pos == "TRUE"

  base_type[is_mhcii & is_ccr2p]  <- "Ccr2+MHCII-hi Mac"
  base_type[is_mhcii & !is_ccr2p] <- "Ccr2-MHCII-hi Mac"

  seuobj_macmonodc$MacType1_final <- factor(
    base_type,
    levels = c(
      "Lyve1+ Mac","Ccr2-MHCII-hi Mac","Ccr2+MHCII-hi Mac",
      "Ccl8+ Mac","Spp1+ Mac","IFN Mac","Fn1+ Mac","Saa3+ Mac","Prolif Mac",
      "Ly6c-hi Mono","Ly6c-lo Mono",
      "cDC1","cDC2","pDC","Migratory DC","Other"
    )
  )
} else {
  seuobj_macmonodc$MacType1_final <- seuobj_macmonodc$MacType1_refined
}

markers_MacType1_final <- find_all_markers_wrapper(
  seu   = seuobj_macmonodc,
  ident = "MacType1_final",
  min_pct = 0.25
)
write.csv(
  markers_MacType1_final,
  file      = file.path("r_objects", "macmonodc_markers_by_MacType1_final.csv"),
  row.names = FALSE
)

save_plot(
  DimPlot(seuobj_macmonodc, reduction="umap",
          group.by="MacType1_final", label=TRUE) +
    ggtitle("MacType1_final") + plot_theme_common(),
  filename = file.path("output_plots", "macmonodc_umap_MacType1_final.pdf"),
  width  = 11,
  height = 7
)


## ---------------------------------------------------------------------------
## 10. Save final Seurat object
## ---------------------------------------------------------------------------

seuobj_macmonodc_final <- seuobj_macmonodc
saveRDS(
  seuobj_macmonodc_final,
  file = file.path("r_objects", "seuobj_macmonodc_final.Rds")
)

message("CD45+ Mac/Mono/DC subcluster pipeline finished. Final object: r_objects/seuobj_macmonodc_final.Rds")

###############################################################################
# End of cd45_MacMonoDc_subcluster.R
###############################################################################
