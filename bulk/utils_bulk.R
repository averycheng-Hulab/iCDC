############################################################
# utils_bulk.R
# Shared utilities for bulk RNA-seq DEG enrichment analysis
# - Gene Ontology (Biological Process) enrichment
# - KEGG pathway enrichment
# - Basic visualization helpers
#
# Species: Mus musculus (org.Mm.eg.db)
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(ggplot2)
  library(enrichplot)
})

############################################################
# 1. Core enrichment function: GO + KEGG
############################################################

enrich_GO_KEGG_mouse <- function(
  genes,
  ntop           = 3000,
  go_ont         = "BP",
  go_padj_cutoff = 0.05,
  go_q_cutoff    = 0.05,
  kegg_p_cutoff  = 0.10,
  kegg_q_cutoff  = 0.10,
  orgdb          = org.Mm.eg.db
) {
  # Ensure character and unique
  genes <- unique(as.character(genes))
  if (length(genes) == 0) {
    return(list(GO = data.frame(), KEGG = data.frame()))
  }

  # Limit to top ntop genes (typically already ordered by padj or p-value)
  genes <- genes[seq_len(min(length(genes), ntop))]

  ## ---- GO enrichment (Biological Process) ----
  ego <- tryCatch(
    enrichGO(
      gene          = genes,
      OrgDb         = orgdb,
      keyType       = "SYMBOL",
      ont           = go_ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = go_padj_cutoff,
      qvalueCutoff  = go_q_cutoff
    ),
    error = function(e) NULL
  )

  if (!is.null(ego) && nrow(ego@result) > 0) {
    ego_simp <- tryCatch(
      simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min),
      error = function(e) ego  # fall back to raw result if simplify fails
    )
    go_res <- as.data.frame(ego_simp@result)
  } else {
    go_res <- data.frame()
  }

  ## ---- KEGG enrichment ----
  genes_id <- tryCatch(
    bitr(
      genes,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = orgdb
    ),
    error = function(e) data.frame()
  )

  if (nrow(genes_id) > 0) {
    ekegg <- tryCatch(
      enrichKEGG(
        gene          = genes_id$ENTREZID,
        keyType       = "kegg",
        organism      = "mmu",
        pvalueCutoff  = kegg_p_cutoff,
        pAdjustMethod = "BH",
        qvalueCutoff  = kegg_q_cutoff
      ),
      error = function(e) NULL
    )

    if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
      kegg_res <- as.data.frame(ekegg@result)

      # Convert ENTREZ IDs back to SYMBOL for readability
      kegg_res$symbol_vec <- vapply(
        strsplit(kegg_res$geneID, "/"),
        FUN.VALUE = character(1),
        FUN = function(x) {
          sym <- tryCatch(
            bitr(x, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = orgdb),
            error = function(e) data.frame()
          )
          if (nrow(sym) == 0) {
            return(NA_character_)
          }
          paste(sym$SYMBOL, collapse = "/")
        }
      )

      # Clean species suffix in pathway description
      kegg_res$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_res$Description)
    } else {
      kegg_res <- data.frame()
    }
  } else {
    kegg_res <- data.frame()
  }

  return(list(GO = go_res, KEGG = kegg_res))
}

############################################################
# 2. Visualization helpers
############################################################

# Barplot for top GO terms (Biological Process)
plot_GO_bar <- function(go_df, n = 10, title = "Top GO terms") {
  if (is.null(go_df) || nrow(go_df) == 0) {
    warning("GO results are empty; cannot plot.")
    return(NULL)
  }

  df <- go_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = n)

  p <- ggplot(df, aes(x = reorder(Description, Count), y = Count)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = title, x = "GO term (BP)", y = "Gene count") +
    theme_bw()
  return(p)
}

# Dotplot for top KEGG pathways
plot_KEGG_dot <- function(kegg_df, n = 10, title = "Top KEGG pathways") {
  if (is.null(kegg_df) || nrow(kegg_df) == 0) {
    warning("KEGG results are empty; cannot plot.")
    return(NULL)
  }

  df <- kegg_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = n)

  p <- ggplot(df, aes(x = Count, y = reorder(Description, Count), size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = title, x = "Gene count", y = "KEGG pathway") +
    theme_bw()
  return(p)
}
