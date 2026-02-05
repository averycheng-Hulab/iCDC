###############################################################################
# utils_scRNA.R
#
# Shared utility functions for all scRNA-seq pipelines:
# - whole_heart_main.R
# - cd45_main.R
# - cd4_T_main.R
# - fibroblast / MacMonoDC / TNK / Neutrophil / B-cell subclustering
# - enrichment modules (indirectly)
#
# Provides standardized:
#   - QC wrappers
#   - SCTransform → PCA → Harmony pipeline
#   - Clustering pipeline
#   - Marker detection functions
#   - Annotation mapping helpers
#   - Plotting helpers
#
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)


###############################################################################
# 1. Common Plot Theme & Save Function
###############################################################################

plot_theme_common <- function() {
  theme_bw() +
    theme(
      axis.title.x = element_text(color="black", face="bold", size=14),
      axis.text.x  = element_text(color="black", face="bold", size=10),
      axis.text.y  = element_text(color="black", face="bold", size=10),
      axis.title.y = element_text(color="black", face="bold", size=14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color="black", size=1),
      plot.title = element_text(hjust=0.5, face="bold", size=16)
    )
}

save_plot <- function(p, filename, width=7, height=5) {
  ggsave(filename, p, width = width, height = height)
}


###############################################################################
# 2. SCTransform → PCA → Harmony Wrapper
###############################################################################
# Used in:
#   - whole_heart_main
#   - cd45_main
#   - cd4_T_main
#   - ALL subclustering pipelines
###############################################################################

run_sct_pca_harmony <- function(
  seu,
  vars_to_regress = "CC.Difference",
  npcs = 50,
  dims_use = 1:40,
  harmony_group = "orig.ident",
  max_iter = 20
){
  message("Running SCTransform...")
  seu <- SCTransform(seu, vars.to.regress = vars_to_regress, verbose = FALSE)

  message("Running PCA...")
  seu <- RunPCA(
    seu, 
    features = VariableFeatures(seu),
    npcs = npcs
  )

  message("Running Harmony...")
  seu <- RunHarmony(
    seu,
    group.by.vars = harmony_group,
    assay.use = "SCT",
    plot_convergence = TRUE,
    max.iter.harmony = max_iter
  )

  return(seu)
}


###############################################################################
# 3. Clustering Wrapper (Neighbors → Clusters → UMAP/TSNE)
###############################################################################

run_clustering <- function(
  seu,
  dims_use = 1:40,
  resolution = 1.0
){
  message("Finding neighbors & clusters...")
  seu <- FindNeighbors(seu, reduction="harmony", dims=dims_use)
  seu <- FindClusters(seu, reduction="harmony", resolution = resolution, dims=dims_use)

  message("Running UMAP/TSNE...")
  seu <- RunTSNE(seu, reduction="harmony", dims=dims_use)
  seu <- RunUMAP(seu, reduction="harmony", dims=dims_use)

  return(seu)
}


###############################################################################
# 4. Marker Detection Wrapper
###############################################################################

find_all_markers_wrapper <- function(
  seu,
  ident = "seurat_clusters",
  min_pct = 0.25
){
  Idents(seu) <- ident
  markers <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = min_pct
  )
  return(markers)
}


###############################################################################
# 5. Annotation Utilities
###############################################################################
# Used in:
#   - Fibroblast Type mapping
#   - MacMonoDC supervised annotation (MacType1 / refined / final)
#   - TNK / Neutrophil / B-cell mapping
###############################################################################

annotate_by_cluster <- function(seu, mapping, new_col) {
  clust_chr <- as.character(seu$seurat_clusters)
  names(mapping) <- as.character(names(mapping))

  mapped <- vapply(
    clust_chr,
    FUN.VALUE = character(1),
    FUN = function(x){
      if (x %in% names(mapping)) mapping[[x]] else "Other"
    }
  )

  levs <- unique(c(unname(mapping), "Other"))
  seu[[new_col]] <- factor(mapped, levels = levs)
  return(seu)
}

# Merge refined annotation back to parent object
merge_annotation_back <- function(
  parent_obj,
  refined_obj,
  parent_col = "MacType1",
  refined_col = "MacType_sel",
  new_col = "MacType1_refined"
){
  meta_main <- parent_obj@meta.data %>% mutate(CellId = rownames(.))
  meta_sel  <- refined_obj@meta.data %>% mutate(CellId = rownames(.))

  merge_ref <- meta_main %>%
    dplyr::select(CellId, !!parent_col) %>%
    left_join(meta_sel %>% dplyr::select(CellId, !!refined_col), by = "CellId") %>%
    mutate(
      !!new_col := if_else(!is.na(.data[[refined_col]]),
                           .data[[refined_col]],
                           .data[[parent_col]])
    )

  rownames(merge_ref) <- merge_ref$CellId
  parent_obj[[new_col]] <- factor(merge_ref[[new_col]])

  return(parent_obj)
}


###############################################################################
# 6. Cell cycle scoring wrapper
###############################################################################
# Given a QC'ed Seurat object (after removing mito/ribo genes),
# and s_genes / g2m_genes (character vectors of gene symbols),
# add S.Score, G2M.Score, Phase, CC.Difference to meta.data and return.
###############################################################################

compute_cell_cycle_scores <- function(
  seu,
  s_genes,
  g2m_genes
){
  message("Computing cell cycle scores...")

  # Temporary object for scoring (SCTransform + PCA)
  tmp <- SCTransform(seu, verbose = FALSE)
  tmp <- RunPCA(tmp, features = VariableFeatures(tmp), verbose = FALSE)
  tmp <- CellCycleScoring(
    tmp,
    s.features   = s_genes,
    g2m.features = g2m_genes,
    set.ident    = TRUE
  )

  # Write scores back to the input object
  seu$S.Score   <- tmp$S.Score
  seu$G2M.Score <- tmp$G2M.Score
  seu$Phase     <- tmp$Phase
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score

  return(seu)
}


###############################################################################
# 7. Reusable Filtering Helpers
###############################################################################

filter_mito_ribo <- function(
  seu,
  mito_pattern = "^mt-",
  ribo_pattern = "^Rp[sl][0-9]"
) {
  counts <- GetAssayData(seu, slot = "counts")
  keep <- !grepl(ribo_pattern, rownames(counts)) &
          !grepl(mito_pattern, rownames(counts))
  seu <- subset(seu, features = rownames(counts)[keep])
  return(seu)
}


###############################################################################
# 8. Per-sample QC + DoubletFinder (shared by all main pipelines)
###############################################################################

process_sample_scRNA <- function(
  matrix_dir,
  sample_name,
  treatment_label,
  min_features = 300,
  min_counts   = 500,
  max_mito     = 5,
  hb_thr       = 5,
  pc_dims      = 1:40,
  doublet_rate = 0.008
){
  message("Processing sample: ", sample_name)

  # 1) Read 10X + basic object
  x <- Read10X(matrix_dir)
  seuobj <- CreateSeuratObject(
    counts       = x,
    project      = sample_name,
    min.cells    = 3,
    min.features = 200
  )
  seuobj$orig.ident <- sample_name
  seuobj$Treatment  <- treatment_label

  # 2) QC metrics
  seuobj[["mito.per"]] <- PercentageFeatureSet(seuobj, pattern = "^mt-|^MT-")
  seuobj[["ribo.per"]] <- PercentageFeatureSet(seuobj, pattern = "^Rp[sl][[:digit:]]")
  seuobj[["hb.per"]]   <- PercentageFeatureSet(seuobj, pattern = "^Hb[ba]-|^HB[AB]")

  # Violin before filtering
  v_pre <- VlnPlot(
    seuobj,
    features = c("nFeature_RNA", "nCount_RNA", "mito.per", "hb.per"),
    pt.size  = 0.01,
    ncol     = 4
  ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  ggsave(
    filename = file.path("QC", paste0("QC_filter_", sample_name, "_before.pdf")),
    plot     = v_pre,
    width    = 7,
    height   = 5
  )

  # 3) Dynamic upper cutoffs (top 1% features, top 2% counts)
  feat_sorted <- sort(seuobj@meta.data$nFeature_RNA, decreasing = TRUE)
  maxFeature  <- feat_sorted[max(1, round(length(feat_sorted) * 0.01))]
  cnt_sorted  <- sort(seuobj@meta.data$nCount_RNA, decreasing = TRUE)
  maxCount    <- cnt_sorted[max(1, round(length(cnt_sorted) * 0.02))]

  seuobj_filt <- subset(
    seuobj,
    subset =
      nFeature_RNA >= min_features &
      nFeature_RNA <  maxFeature   &
      nCount_RNA    >= min_counts  &
      nCount_RNA    <  maxCount    &
      hb.per        <  hb_thr      &
      mito.per      <= max_mito
  )

  # Violin after filtering
  v_post <- VlnPlot(
    seuobj_filt,
    features = c("nFeature_RNA", "nCount_RNA", "mito.per", "hb.per"),
    pt.size  = 0.01,
    ncol     = 4
  ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  ggsave(
    filename = file.path("QC", paste0("QC_filter_", sample_name, "_after.pdf")),
    plot     = v_post,
    width    = 7,
    height   = 5
  )

  # 4) Prepare for DoubletFinder
  seuobj_df <- NormalizeData(seuobj_filt, verbose = FALSE) |>
    FindVariableFeatures(nfeatures = 3000, verbose = FALSE) |>
    ScaleData(verbose = FALSE)

  seuobj_df <- RunPCA(seuobj_df, verbose = FALSE)
  seuobj_df <- RunTSNE(seuobj_df, dims = pc_dims, verbose = FALSE)
  seuobj_df <- RunUMAP(seuobj_df, dims = pc_dims, verbose = FALSE)
  seuobj_df <- FindNeighbors(seuobj_df, dims = pc_dims, verbose = FALSE)
  seuobj_df <- FindClusters(seuobj_df, verbose = FALSE)

  # 5) DoubletFinder parameter sweep (v2 / v3 自动兼容)
  if ("paramSweep_v3" %in% ls(getNamespace("DoubletFinder"))) {
    sweep_res  <- DoubletFinder::paramSweep_v3(seuobj_df, PCs = pc_dims, sct = FALSE)
    sweep_stat <- DoubletFinder::summarizeSweep(sweep_res, GT = FALSE)
    bcmvn      <- DoubletFinder::find.pK(sweep_stat)
  } else {
    sweep_res  <- DoubletFinder::paramSweep(seuobj_df, PCs = pc_dims, sct = FALSE)
    sweep_stat <- DoubletFinder::summarizeSweep(sweep_res, GT = FALSE)
    bcmvn      <- DoubletFinder::find.pK(sweep_stat)
  }

  pK_best <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  homotypic_prop <- DoubletFinder::modelHomotypic(seuobj_df$seurat_clusters)
  DoubleRate     <- (ncol(seuobj_df) / 1000) * doublet_rate
  nExp_poi       <- round(DoubleRate * nrow(seuobj_df@meta.data))
  nExp_poi_adj   <- round(nExp_poi * (1 - homotypic_prop))

  if ("doubletFinder_v3" %in% ls(getNamespace("DoubletFinder"))) {
    seuobj_df <- DoubletFinder::doubletFinder_v3(
      seuobj_df,
      PCs        = pc_dims,
      pN         = 0.25,
      pK         = pK_best,
      nExp       = nExp_poi_adj,
      reuse.pANN = FALSE,
      sct        = FALSE
    )
    df_col <- grep("DF.classifications", colnames(seuobj_df@meta.data), value = TRUE)
  } else {
    seuobj_df <- DoubletFinder::doubletFinder(
      seuobj_df,
      PCs        = pc_dims,
      pN         = 0.25,
      pK         = pK_best,
      nExp       = nExp_poi_adj,
      reuse.pANN = FALSE,
      sct        = FALSE
    )
    df_col <- grep("doublet|Doublet|DF.class", colnames(seuobj_df@meta.data), value = TRUE)[1]
  }

  # 6) Doublet UMAP plot
  if (length(df_col) > 0 && !is.na(df_col[1]) && nzchar(df_col[1])) {
    p_doublet <- DimPlot(seuobj_df, group.by = df_col[1], reduction = "umap") +
      ggtitle(paste0(sample_name, " doublet classification"))
    ggsave(
      filename = file.path("QC", paste0("QC_doublet_", sample_name, ".pdf")),
      plot     = p_doublet,
      width    = 8,
      height   = 4
    )
    seuobj_sing <- subset(seuobj_df, subset = get(df_col[1]) == "Singlet")
  } else {
    warning("No DoubletFinder classification column found for sample ", sample_name, ". Using all cells.")
    seuobj_sing <- seuobj_df
  }

  seuobj_sing$orig.ident <- sample_name
  seuobj_sing$Treatment  <- treatment_label

  return(seuobj_sing)
}


###############################################################################
# End of utils_scRNA.R
###############################################################################
