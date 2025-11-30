###############################################################################
# whole_heart_fibroblast_ecm.R
#
# Fibroblast ECM Regulator gene program analysis (only ECM_Regulator_Score)
#
# - Input:
#     r_objects/seuobj_fibroblast.Rds
#       (output from whole_heart_fibroblast.R)
#
# - Steps:
#     1) Load fibroblast Seurat object
#     2) Load ECM_regulators.tsv
#     3) Compute ECM_Regulator_Score via AddModuleScore
#     4) Visualize ECM_Regulator_Score on UMAP
#     5) Violin plots of ECM_Regulator_Score by FibType and Treatment
#
# - Output:
#     r_objects/seuobj_fibroblast_ecm.Rds
###############################################################################

## ---------------------------------------------------------------------------
## 0. Packages, utils, options
## ---------------------------------------------------------------------------

required_pkgs <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "stringr"
)

missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Please install them before running this script."
  )
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

source("utils_scRNA.R")

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)

dir.create("fibroblast_plots", showWarnings = FALSE)
dir.create("r_objects", showWarnings = FALSE)


## ---------------------------------------------------------------------------
## 1. Load fibroblast object
## ---------------------------------------------------------------------------

seuobj_fibro <- readRDS(file.path("r_objects", "seuobj_fibroblast.Rds"))

if (!"FibType" %in% colnames(seuobj_fibro@meta.data)) {
  stop("FibType column not found. Run whole_heart_fibroblast.R first.")
}

if ("SCT" %in% names(seuobj_fibro@assays)) {
  DefaultAssay(seuobj_fibro) <- "SCT"
} else {
  DefaultAssay(seuobj_fibro) <- "RNA"
}


## ---------------------------------------------------------------------------
## 2. Load ECM Regulator list only
## ---------------------------------------------------------------------------

reg_file <- file.path("ECM_genelist", "ECMregulators.tsv")

if (!file.exists(reg_file)) {
  stop("ECMregulators.tsv not found in ECM_genelist/.")
}

reg_frame <- read.table(
  reg_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

reg_genes_raw <- reg_frame[[1]]
reg_genes_std <- unique(str_to_title(reg_genes_raw))
reg_genes_in_obj <- intersect(reg_genes_std, rownames(seuobj_fibro))

if (length(reg_genes_in_obj) == 0) {
  stop("None of ECM regulator genes found in Seurat object.")
}

message(length(reg_genes_in_obj), " ECM regulator genes detected.")


## ---------------------------------------------------------------------------
## 3. Compute ECM_Regulator_Score (single module score)
## ---------------------------------------------------------------------------

seuobj_fibro <- AddModuleScore(
  seuobj_fibro,
  features = list(reg_genes_in_obj),
  name     = "ECM_Regulator"
)

# AddModuleScore names new column like 'ECM_Regulator1'
new_col <- tail(colnames(seuobj_fibro@meta.data), 1)
colnames(seuobj_fibro@meta.data)[colnames(seuobj_fibro@meta.data) == new_col] <- "ECM_Regulator_Score"


## ---------------------------------------------------------------------------
## 4. UMAP visualization
## ---------------------------------------------------------------------------

if ("umap" %in% names(seuobj_fibro@reductions)) {
  p_umap <- FeaturePlot(
    seuobj_fibro,
    features  = "ECM_Regulator_Score",
    reduction = "umap"
  ) +
    ggtitle("Fibroblasts: ECM_Regulator_Score (UMAP)") +
    plot_theme_common()

  save_plot(
    p_umap,
    file.path("fibroblast_plots", "fibro_ECM_Regulator_UMAP.pdf"),
    width = 8,
    height = 6
  )
} else {
  warning("UMAP not found. Skipping FeaturePlot.")
}


## ---------------------------------------------------------------------------
## 5. Violin plots
## ---------------------------------------------------------------------------

# 5.1 By FibType
p_fib <- VlnPlot(
  seuobj_fibro,
  features = "ECM_Regulator_Score",
  group.by = "FibType",
  pt.size  = 0
) +
  ggtitle("ECM_Regulator_Score by FibType") +
  plot_theme_common()

save_plot(
  p_fib,
  file.path("fibroblast_plots", "fibro_ECM_Regulator_by_FibType.pdf"),
  width = 9,
  height = 6
)

# 5.2 By Treatment if available
if ("Treatment" %in% colnames(seuobj_fibro@meta.data)) {
  p_treat <- VlnPlot(
    seuobj_fibro,
    features = "ECM_Regulator_Score",
    group.by = "Treatment",
    pt.size  = 0
  ) +
    ggtitle("ECM_Regulator_Score by Treatment") +
    plot_theme_common()

  save_plot(
    p_treat,
    file.path("fibroblast_plots", "fibro_ECM_Regulator_by_Treatment.pdf"),
    width = 7,
    height = 6
  )
}


## ---------------------------------------------------------------------------
## 6. Save updated object
## ---------------------------------------------------------------------------

saveRDS(
  seuobj_fibro,
  file.path("r_objects", "seuobj_fibroblast_ecm.Rds")
)

message("ECM regulator analysis finished. Saved to r_objects/seuobj_fibroblast_ecm.Rds")

###############################################################################
# End of whole_heart_fibroblast_ecm.R
###############################################################################
