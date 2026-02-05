###############################################################################
# enrich_scRNA.R
#
# Unified DEG + GO/KEGG enrichment workflow for any scRNA-seq cell type.
#
# Usage:
#   source("enrich_scRNA.R")
#   run_enrich_scRNA(
#        obj      = seu_object,
#        celltype = "F-SL",
#        outdir   = "Fib_SL_enrich"
#   )
#
# Functionality:
#   1) IR vs Sham DEG (SCT)
#   2) IR vs DC DEG (SCT)
#   3) Intersection of up/down DEGs
#   4) GO/KEGG enrichment
#   5) Save DEG + enrichment tables
#   6) Bubble plots for GO and KEGG (Up / Down)
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(writexl)
})

# Load utilities: theme + save_plot
source("utils_scRNA.R")

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)

###############################################################################
# 1. Differential Expression: IR vs Sham, IR vs DC
###############################################################################

get_two_comparisons <- function(obj, celltype) {
  Idents(obj) <- "Treatment"

  deg1 <- FindMarkers(
    obj,
    ident.1         = "IR",
    ident.2         = "Sham",
    assay           = "SCT",
    subset.ident    = celltype,
    logfc.threshold = 0.3,
    min.pct         = 0.5
  )
  deg1$gene <- rownames(deg1)

  deg2 <- FindMarkers(
    obj,
    ident.1         = "IR",
    ident.2         = "DC",
    assay           = "SCT",
    subset.ident    = celltype,
    logfc.threshold = 0.3,
    min.pct         = 0.5
  )
  deg2$gene <- rownames(deg2)

  list(IR_vs_Sham = deg1, IR_vs_DC = deg2)
}

###############################################################################
# 2. Extract shared up/down DEGs
###############################################################################

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

###############################################################################
# 3. GO / KEGG enrichment wrapper
###############################################################################

run_GO_KEGG <- function(gene_vec, top = 100) {
  gene_vec <- gene_vec[1:min(top, length(gene_vec))]
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

###############################################################################
# 4. Bubble plot (same style as cd4 enrichment)
###############################################################################

plot_bubble <- function(df, title_text, file_out) {
  if (is.null(df) || nrow(df) == 0) {
    return(invisible(NULL))
  }

  df_plot <- df[1:min(10, nrow(df)), ]
  df_plot$Description <- factor(df_plot$Description,
                                levels = rev(df_plot$Description))

  p <- ggplot(df_plot, aes(x = 1, y = Description)) +
    geom_point(aes(size = Count, color = -log10(pvalue))) +
    scale_color_gradient(low = "#56B1F7", high = "#CA0020") +
    scale_size(range = c(3, 8)) +
    labs(x = "", y = "", title = title_text) +
    plot_theme_common()

  save_plot(p, file_out, width = 6, height = 5)
}

###############################################################################
# 5. Master function: run_enrich_scRNA
###############################################################################

run_enrich_scRNA <- function(
  obj,
  celltype,
  outdir = paste0("enrich_", celltype),
  top = 100
) {
  dir.create(outdir, showWarnings = FALSE)

  comp <- get_two_comparisons(obj, celltype)
  deg_shared <- extract_intersection_deg(
    comp$IR_vs_Sham,
    comp$IR_vs_DC
  )

  enrich_up   <- run_GO_KEGG(deg_shared$up,   top = top)
  enrich_down <- run_GO_KEGG(deg_shared$down, top = top)

  # Save DEG tables
  write_xlsx(
    list(
      IR_vs_Sham  = comp$IR_vs_Sham,
      IR_vs_DC    = comp$IR_vs_DC,
      shared_up   = data.frame(gene = deg_shared$up),
      shared_down = data.frame(gene = deg_shared$down)
    ),
    file.path(outdir, "DEG_tables.xlsx")
  )

  # Save GO/KEGG tables
  write_xlsx(
    list(
      GO_up    = enrich_up$GO,
      KEGG_up  = enrich_up$KEGG,
      GO_down  = enrich_down$GO,
      KEGG_down= enrich_down$KEGG
    ),
    file.path(outdir, "GO_KEGG_tables.xlsx")
  )

  # Bubble plots
  if (!is.null(enrich_up$GO))
    plot_bubble(
      enrich_up$GO,
      paste0(celltype, ": GO Up"),
      file.path(outdir, "GO_up_bubble.pdf")
    )

  if (!is.null(enrich_down$GO))
    plot_bubble(
      enrich_down$GO,
      paste0(celltype, ": GO Down"),
      file.path(outdir, "GO_down_bubble.pdf")
    )

  if (!is.null(enrich_up$KEGG))
    plot_bubble(
      enrich_up$KEGG,
      paste0(celltype, ": KEGG Up"),
      file.path(outdir, "KEGG_up_bubble.pdf")
    )

  if (!is.null(enrich_down$KEGG))
    plot_bubble(
      enrich_down$KEGG,
      paste0(celltype, ": KEGG Down"),
      file.path(outdir, "KEGG_down_bubble.pdf")
    )

  invisible(
    list(
      DEG   = deg_shared,
      GO    = list(up = enrich_up$GO,   down = enrich_down$GO),
      KEGG  = list(up = enrich_up$KEGG, down = enrich_down$KEGG)
    )
  )
}

###############################################################################
# End of enrich_scRNA.R
###############################################################################
