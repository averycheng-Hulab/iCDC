############################################################
# bulk_1-2.R
# Bulk RNA-seq DEG enrichment (T cell / DC bulk)
#
# - Input: Excel DEG tables (FPKM-based)
#          columns must contain at least:
#          - gene_symbol
#          - log2FoldChange
#          - padj
#          - two expression columns used for exp_cutoff 
#
# - Output: DEG up/down lists, GO/KEGG enrichment tables,
#           basic example plots
#
# - Datasets:
#   bulk_1: T cell bulk RNA-seq
#   bulk_2: DC bulk RNA-seq
#
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# Source shared utilities
source(file.path("bulk","utils","utils_bulk.R"))

############################################################
# 1. Load DEG tables
############################################################

# Input directory containing DEG Excel files (FPKM-based)
input_dir  <- file.path("bulk", "input", "bulk_1-2")
output_dir <- file.path("bulk", "output", "bulk_1-2")
plot_dir   <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

deg_files <- c(
  file.path(input_dir, "DEG_Tcell_example.xlsx"),
  file.path(input_dir, "DEG_DC_example.xlsx")
)

# Names corresponding to each comparison
deg_names <- c("Tcell_Condition_vs_Control", "DC_Condition_vs_Control")

# Check consistency
stopifnot(length(deg_files) == length(deg_names))

# Parameters
exp_cutoff  <- 35   # minimal expression cutoff (applied to 2 average expression columns of two groups)
logfc_cutoff <- 1   # absolute log2FC threshold
ntop        <- 3000 # maximum number of genes for enrichment

# e.g. column x: average value of exp_group1; column y: average value of exp_group2
expr_col_idx <- c(x, y)

DEG_list <- lapply(deg_files, function(f) {
  message("Reading: ", f)
  as.data.frame(read_xlsx(f))
})
names(DEG_list) <- deg_names

############################################################
# 2. GO/KEGG analysis and results display
############################################################

DEG_up_list       <- list()
DEG_down_list     <- list()
DEG_up_GO_list    <- list()
DEG_down_GO_list  <- list()
DEG_up_KEGG_list  <- list()
DEG_down_KEGG_list<- list()

for (i in seq_along(DEG_list)) {
  comp_name <- deg_names[i]
  df <- DEG_list[[i]]

  message("Processing comparison: ", comp_name)

  # Sanity checks
  required_cols <- c("gene_symbol", "log2FoldChange", "padj")
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns in ", comp_name, ": ",
         paste(setdiff(required_cols, colnames(df)), collapse = ", "))
  }

  # Expression filter using two columns (as in original FPKM-based code)
  # We assume expr_col_idx points to expression columns
  idx_exp <- which(df[, expr_col_idx[1]] > exp_cutoff |
                   df[, expr_col_idx[2]] > exp_cutoff)
  df <- df[idx_exp, , drop = FALSE]

  # Upregulated genes
  up_genes <- df %>%
    dplyr::filter(log2FoldChange > logfc_cutoff) %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene_symbol) %>%
    unique()

  # Downregulated genes
  down_genes <- df %>%
    dplyr::filter(log2FoldChange < -logfc_cutoff) %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene_symbol) %>%
    unique()

  DEG_up_list[[comp_name]]   <- up_genes
  DEG_down_list[[comp_name]] <- down_genes

  # Enrichment analysis
  up_enrich   <- enrich_GO_KEGG_mouse(up_genes, ntop = ntop)
  down_enrich <- enrich_GO_KEGG_mouse(down_genes, ntop = ntop)

  DEG_up_GO_list[[comp_name]]    <- up_enrich$GO
  DEG_up_KEGG_list[[comp_name]]  <- up_enrich$KEGG
  DEG_down_GO_list[[comp_name]]  <- down_enrich$GO
  DEG_down_KEGG_list[[comp_name]]<- down_enrich$KEGG


  # plots (upregulated GO/KEGG)
  p_go_bar <- plot_enrich_bar_neglogp(
    up_enrich$GO,
    n = 10,
    title = paste0("Top10 GO (Up) - ", comp_name),
    p_col = "pvalue"   
  )

  p_kegg_bar <- plot_enrich_bar_neglogp(
    up_enrich$KEGG,
    n = 10,
    title = paste0("Top10 KEGG (Up) - ", comp_name),
    p_col = "pvalue"   # or "p.adjust"
  )

  if (!is.null(p_go_bar)) {
    ggsave(file.path(plot_dir, paste0("GO_up_", comp_name, ".pdf")),
           p_go_bar, width = 6, height = 4)
  }
  if (!is.null(p_kegg_bar)) {
    ggsave(file.path(plot_dir, paste0("KEGG_up_", comp_name, ".pdf")),
           p_kegg_bar, width = 6, height = 4)
  }

  #####combine
  if (!is.null(p_go_bar) || !is.null(p_kegg_bar)) {
    p_combo <- patchwork::wrap_plots(
      list(p_go_bar, p_kegg_bar),
      ncol = 1
    ) +
      patchwork::plot_annotation(
        title = paste0("Enrichment (up) - ", comp_name)
      )

    ggsave(
      filename = file.path(plot_dir, paste0("Enrich_up_top10_GO_KEGG_", comp_name, ".pdf")),
      plot = p_combo, width = 6, height = 9
    )
  }

  # plots (downregulated GO/KEGG)
  p_go_bar <- plot_enrich_bar_neglogp(
    down_enrich$GO,
    n = 10,
    title = paste0("Top10 GO (Down) - ", comp_name),
    p_col = "pvalue"   
  )

  p_kegg_bar <- plot_enrich_bar_neglogp(
    down_enrich$KEGG,
    n = 10,
    title = paste0("Top10 KEGG (Down) - ", comp_name),
    p_col = "pvalue"   # or "p.adjust"
  )

  if (!is.null(p_go_bar)) {
    ggsave(file.path(plot_dir, paste0("GO_down_", comp_name, ".pdf")),
           p_go_bar, width = 6, height = 4)
  }
  if (!is.null(p_kegg_bar)) {
    ggsave(file.path(plot_dir, paste0("KEGG_down_", comp_name, ".pdf")),
           p_kegg_bar, width = 6, height = 4)
  }

  #####combine
  if (!is.null(p_go_bar) || !is.null(p_kegg_bar)) {
    p_combo <- patchwork::wrap_plots(
      list(p_go_bar, p_kegg_bar),
      ncol = 1
    ) +
      patchwork::plot_annotation(
        title = paste0("Enrichment (down) - ", comp_name)
      )

    ggsave(
      filename = file.path(plot_dir, paste0("Enrich_down_top10_GO_KEGG_", comp_name, ".pdf")),
      plot = p_combo, width = 6, height = 9
    )
  }
}

############################################################
# 3. intersections between comparisons
############################################################

if (length(DEG_up_list) >= 2) {
  common_up   <- Reduce(intersect, DEG_up_list)
  common_down <- Reduce(intersect, DEG_down_list)

  DEG_up_list[["Common"]]   <- common_up
  DEG_down_list[["Common"]] <- common_down

  common_up_enrich   <- enrich_GO_KEGG_mouse(common_up,   ntop = ntop)
  common_down_enrich <- enrich_GO_KEGG_mouse(common_down, ntop = ntop)

  DEG_up_GO_list[["Common"]]    <- common_up_enrich$GO
  DEG_up_KEGG_list[["Common"]]  <- common_up_enrich$KEGG
  DEG_down_GO_list[["Common"]]  <- common_down_enrich$GO
  DEG_down_KEGG_list[["Common"]]<- common_down_enrich$KEGG
}

############################################################
# 4. Save results
############################################################

write_xlsx(DEG_up_list,
           path = file.path(output_dir, "DEG_up_gene_lists.xlsx"))
write_xlsx(DEG_down_list,
           path = file.path(output_dir, "DEG_down_gene_lists.xlsx"))
write_xlsx(DEG_up_GO_list,
           path = file.path(output_dir, "DEG_up_GO_results.xlsx"))
write_xlsx(DEG_down_GO_list,
           path = file.path(output_dir, "DEG_down_GO_results.xlsx"))
write_xlsx(DEG_up_KEGG_list,
           path = file.path(output_dir, "DEG_up_KEGG_results.xlsx"))
write_xlsx(DEG_down_KEGG_list,
           path = file.path(output_dir, "DEG_down_KEGG_results.xlsx"))

message("Analysis completed. Results in: ", output_dir)
