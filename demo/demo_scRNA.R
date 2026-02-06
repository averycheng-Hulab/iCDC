if (!dir.exists("demo") || !dir.exists("scRNA")) stop("Please run from repository root.")

suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(readr)
  library(writexl)
  library(dplyr)
  library(stringr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
})

# Source shared utilities
source(file.path("scRNA","utils","utils_scRNA.R"))

input_dir  <- file.path("demo","scRNA_input")
output_dir <- file.path("demo","scRNA_output")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)


## 1. Load object

demo_fib <- readRDS(file.path("demo","scRNA_input","demo_fib.rds"))

## 2. Markers (post-annotation)

Idents(demo_fib) <- "FibType3"
markers_post <- find_all_markers_wrapper(demo_fib, ident="FibType3")
write.csv(markers_post, file.path(output_dir,"fibro_markers_by_FibType3.csv"), row.names=FALSE)


## 2. FibType3 visualization

save_plot(
  CellDimPlot(demo_fib, group.by="FibType3", reduction="UMAP",
              label=TRUE, label_insitu=TRUE, pt.size=0.1) +
    plot_theme_common(),
  file.path(output_dir,"fibro_umap_FibType3.pdf"),
  width=5, height=4
)

save_plot(
  CellDimPlot(demo_fib, group.by="FibType3", reduction="TSNE",
              label=TRUE, label_insitu=TRUE, pt.size=0.1) +
    plot_theme_common(),
  file.path(output_dir,"fibro_tsne_FibType3.pdf"),
  width=5, height=4
)

save_plot(
  CellDimPlot(demo_fib, group.by="FibType3", split.by="Treatment",
              reduction="UMAP", label=TRUE, pt.size=0.08),
  file.path(output_dir,"fibro_umap_FibType3_by_treatment.pdf"),
  width=10, height=10
)

save_plot(
  CellStatPlot(demo_fib, stat.by="FibType3", group.by="Treatment",
               plot_type="trend", label=TRUE) +
    plot_theme_common(),
  file.path(output_dir,"fibro_cellstat_by_treatment.pdf"),
  width=5, height=5
)


## 3. ECMReg score
#Load ECM Regulator list

reg_file <- file.path("scRNA","ECM_genelist", "ECMregulators.csv")

if (!file.exists(reg_file)) {
  stop("ECMregulators.csv not found in ECM_genelist/.")
}

reg_frame <- read.csv(
  reg_file,
  header = TRUE
)

reg_genes_raw <- reg_frame[[1]]
reg_genes_std <- unique(str_to_title(reg_genes_raw))
reg_genes_in_obj <- intersect(reg_genes_std, rownames(demo_fib))

if (length(reg_genes_in_obj) == 0) {
  stop("None of ECM regulator genes found in Seurat object.")
}

message(length(reg_genes_in_obj), " ECM regulator genes detected.")


#Compute ECM_Regulator_Score (single module score)

demo_fib <- AddModuleScore(
  demo_fib,
  features = list(reg_genes_in_obj),
  name     = "ECM_Regulator"
)

# AddModuleScore names new column like 'ECM_Regulator1'
new_col <- tail(colnames(demo_fib@meta.data), 1)
colnames(demo_fib@meta.data)[colnames(demo_fib@meta.data) == new_col] <- "ECM_Regulator_Score"


# UMAP visualization1

save_plot(
  FeatureStatPlot(
  srt = demo_fib, group.by = "FibType3", bg.by = "FibType3",
  stat.by = c("ECM_Regulator_Score"), 
  add_box = TRUE,
  comparisons = list(
  c("F-Myo","F-SL"),
  c("F-Myo","F-SH"),
  c("F-Myo","F-Act"),
  c("F-Myo","F-IFNs"),
  c("F-Myo","F-IR")
  ),
  ncol = 1
  #legend.position = "bottom",
)+
   plot_theme_common(),
  file.path(output_dir,"fibro_ECMreg_score_byFibType3.pdf"),
  width=6, height=4
)


# UMAP visualization1

save_plot(
  FeatureStatPlot(
  srt = demo_fib, group.by = "Treatment", bg.by = "Treatment",
  stat.by = c("ECM_Regulator_Score"), 
  add_box = TRUE,
  comparisons = list(
  c("Sham","IR"),
  c("IR","Vector"),
  c("IR","DC"),
  c("Vector","DC")
  ),
  ncol = 1
  #legend.position = "bottom",
)+
   plot_theme_common(),
  file.path(output_dir,"fibro_ECMreg_score_byTreats.pdf"),
  width=6, height=4
)


## 4. GO/KEGG analysis

get_two_comparisons <- function(obj, celltype, celltype_col = "CellType1") {
  keep <- obj[[celltype_col]][, 1] == celltype
  obj_ct <- obj[, keep]

  Idents(obj_ct) <- "Treatment"

  deg1 <- FindMarkers(obj_ct, ident.1="IR", ident.2="Sham",
                      assay="SCT", logfc.threshold=0.3, min.pct=0.5)
  deg1$gene <- rownames(deg1)

  deg2 <- FindMarkers(obj_ct, ident.1="IR", ident.2="DC",
                      assay="SCT", logfc.threshold=0.3, min.pct=0.5)
  deg2$gene <- rownames(deg2)

  list(IR_vs_Sham = deg1, IR_vs_DC = deg2)
}


extract_intersection_deg <- function(deg1, deg2, p_cut = 0.01) {
  up1   <- deg1 %>% filter(avg_log2FC > 0, p_val < p_cut) %>% pull(gene)
  up2   <- deg2 %>% filter(avg_log2FC > 0, p_val < p_cut) %>% pull(gene)

  down1 <- deg1 %>% filter(avg_log2FC < 0, p_val < p_cut) %>% pull(gene)
  down2 <- deg2 %>% filter(avg_log2FC < 0, p_val < p_cut) %>% pull(gene)

  list(
    up   = intersect(up1, up2),
    down = intersect(down1, down2)
  )
}


run_GO_KEGG <- function(gene_vec) {
  gene_vec <- gene_vec[!is.na(gene_vec)]

  if (length(gene_vec) < 1) {
    return(list(GO = NULL, KEGG = NULL))
  }

  # GO enrichment
  ego <- enrichGO(
    gene          = gene_vec,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    GO_df <- simplify(
      ego,
      cutoff     = 0.7,
      by         = "p.adjust",
      select_fun = min
    ) |> as.data.frame()
  } else {
    GO_df <- NULL
  }

  # Convert SYMBOL → ENTREZ
  gene_df <- bitr(
    gene_vec,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  )

  # KEGG enrichment
  if (!is.null(gene_df) && nrow(gene_df) > 0) {
    ekegg <- enrichKEGG(
      gene          = gene_df$ENTREZID,
      organism      = "mmu",
      keyType       = "kegg",
      pvalueCutoff  = 1,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1
    )

    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      kegg_df <- as.data.frame(ekegg@result)

      # Convert KEGG gene IDs back to SYMBOL
      symbol_vec <- character(nrow(kegg_df))
      for (i in seq_len(nrow(kegg_df))) {
        ids <- unlist(strsplit(kegg_df$geneID[i], "/"))
        conv <- bitr(
          ids,
          fromType = "ENTREZID",
          toType   = "SYMBOL",
          OrgDb    = org.Mm.eg.db
        )
        symbol_vec[i] <- paste(conv$SYMBOL, collapse = "/")
      }
      kegg_df$symbol_vec <- symbol_vec

      # Clean description
      kegg_df$Description <- gsub(
        " - Mus musculus \\(house mouse\\)",
        "",
        kegg_df$Description
      )

      if ("category" %in% colnames(kegg_df)) {
        kegg_df <- kegg_df[kegg_df$category != "Human Diseases", ]
      }
    } else {
      kegg_df <- NULL
    }
  } else {
    kegg_df <- NULL
  }

  list(GO = GO_df, KEGG = kegg_df)
}

plot_bubble <- function(df, title_text, file_out) {
  if (is.null(df) || nrow(df) == 0) {
    return(invisible(NULL))
  }

  df_plot <- df[1:min(10, nrow(df)), ]
  df_plot$Description <- stringr::str_wrap(df_plot$Description, width = 35)

  df_plot$Description <- factor(df_plot$Description,
                                levels = rev(df_plot$Description))

  p <- ggplot(df_plot, aes(x = 1, y = Description)) +
    geom_point(aes(size = Count, color = -log10(pvalue))) +
    scale_x_continuous(limits = c(1, 1), expand = c(0, 0))+
    scale_color_gradient(low = "#56B1F7", high = "#CA0020") +
    scale_size(range = c(3, 8)) +
    labs(x = "", y = "", title = title_text) +
    plot_theme_common()

  save_plot(p, file_out, width =6, height = 5)
}


  comp <- get_two_comparisons(demo_fib, "Fibroblast")
  deg_shared <- extract_intersection_deg(
    comp$IR_vs_Sham,
    comp$IR_vs_DC
  )

  enrich_up   <- run_GO_KEGG(deg_shared$up)
  enrich_down <- run_GO_KEGG(deg_shared$down)

  # Save DEG tables
  write_xlsx(
    list(
      IR_vs_Sham  = comp$IR_vs_Sham,
      IR_vs_DC    = comp$IR_vs_DC,
      shared_up   = data.frame(gene = deg_shared$up),
      shared_down = data.frame(gene = deg_shared$down)
    ),
    file.path(output_dir, "DEG_tables.xlsx")
  )

  # Save GO/KEGG enrichment results of I/R-induced genes reversed by iCDC
  write_xlsx(
    list(
      GO_up    = enrich_up$GO,
      KEGG_up  = enrich_up$KEGG
    ),
    file.path(output_dir, "GO_KEGG_tables.xlsx")
  )

  # Bubble plots
  if (!is.null(enrich_up$GO))
    plot_bubble(
      enrich_up$GO[complete.cases(enrich_up$GO),],
      paste0("Fibroblat", ": GO Up in IR"),
      file.path(output_dir, "GO_up_bubble.pdf")
    )


  if (!is.null(enrich_up$KEGG))
    plot_bubble(
      enrich_up$KEGG[complete.cases(enrich_up$KEGG),],
      paste0("Fibroblast", ": KEGG Up in IR"),
      file.path(output_dir, "KEGG_up_bubble.pdf")
    )


  invisible(
    list(
      DEG   = deg_shared,
      GO    = list(up = enrich_up$GO),
      KEGG  = list(up = enrich_up$KEGG)
    )
  )

saveRDS(demo_fib, file = file.path(output_dir, "demo_fib2.rds"), compress = "xz" )
message("demo_scRNA analysis completed. Results in: ", output_dir)






























