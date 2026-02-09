if (!dir.exists("demo") || !dir.exists("bulk")) stop("Please run from repository root.")

suppressPackageStartupMessages({
  library(writexl)
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(patchwork)
})

# Source shared utilities
source(file.path("bulk","utils","utils_bulk.R"))

############################################################
# 1. Load DEG tables
############################################################

input_dir  <- file.path("demo","bulk_input")
output_dir <- file.path("demo","bulk_output")
plot_dir   <- output_dir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# DEG files 
deg_files <- c(
  "CAR_vs_IR.csv",  
  "CAR_vs_Vector.csv"
)


contrast_names <- c("CAR_vs_IR", "CAR_vs_Vector")

deg_paths <- file.path(input_dir, deg_files)
names(deg_paths) <- contrast_names

# Filtering parameters for DEGs
logfc_cutoff  <- 1      # |log2FC| > 1
pvalue_cutoff <- 0.01   # p-value threshold
ntop          <- 3000


DEG_list <- lapply(deg_paths, function(f) {
  message("Reading: ", f)
  read.csv(f, header = TRUE, stringsAsFactors = FALSE)
})
names(DEG_list) <- contrast_names


############################################################
# 2. GO/KEGG analysis and results display
############################################################

DEG_up_list        <- list()
DEG_down_list      <- list()
DEG_up_GO_list     <- list()
DEG_down_GO_list   <- list()
DEG_up_KEGG_list   <- list()
DEG_down_KEGG_list <- list()

for (comp_name in names(DEG_list)) {
  df <- DEG_list[[comp_name]]
  message("Processing comparison: ", comp_name)

  # Sanity checks
  required_cols <- c("gene_symbol", "log2FoldChange", "pvalue", "padj")
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns in ", comp_name, ": ",
         paste(setdiff(required_cols, colnames(df)), collapse = ", "))
  }

  # Upregulated DEGs: log2FC > 1 & pvalue < 0.01
  up_df <- df %>%
    dplyr::filter(log2FoldChange >  logfc_cutoff,
                  pvalue         <  pvalue_cutoff) %>%
    dplyr::arrange(padj)

  # Downregulated DEGs: log2FC < -1 & pvalue < 0.01
  down_df <- df %>%
    dplyr::filter(log2FoldChange < -logfc_cutoff,
                  pvalue         <  pvalue_cutoff) %>%
    dplyr::arrange(padj)

  DEG_up_list[[comp_name]]   <- up_df
  DEG_down_list[[comp_name]] <- down_df

  # Extract gene symbols (remove potential spaces)
  up_genes <- gsub(" ", "", up_df$gene_symbol)
  down_genes <- gsub(" ", "", down_df$gene_symbol)

  # Enrichment
  up_enrich   <- enrich_GO_KEGG_mouse(up_genes,   ntop = ntop)
  down_enrich <- enrich_GO_KEGG_mouse(down_genes, ntop = ntop)

  DEG_up_GO_list[[comp_name]]     <- up_enrich$GO
  DEG_up_KEGG_list[[comp_name]]   <- up_enrich$KEGG
  DEG_down_GO_list[[comp_name]]   <- down_enrich$GO
  DEG_down_KEGG_list[[comp_name]] <- down_enrich$KEGG


  # plots (upregulated GO/KEGG)
  p_go_bar <- plot_enrich_bar_neglogp(
    up_enrich$GO,
    n = 10,
    title = paste0("Top10 GO (Up) - ", comp_name),
    p_col = "pvalue"   # or "p.adjust" if you prefer adjusted p
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
    p_col = "pvalue"   # or "p.adjust" if you prefer adjusted p
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
# 3. Save results
############################################################

write_xlsx(DEG_up_list,
           path = file.path(output_dir, "DEG_up_tables.xlsx"))
write_xlsx(DEG_down_list,
           path = file.path(output_dir, "DEG_down_tables.xlsx"))

# GO / KEGG results
write_xlsx(DEG_up_GO_list,
           path = file.path(output_dir, "DEG_up_GO_results.xlsx"))
write_xlsx(DEG_down_GO_list,
           path = file.path(output_dir, "DEG_down_GO_results.xlsx"))
write_xlsx(DEG_up_KEGG_list,
           path = file.path(output_dir, "DEG_up_KEGG_results.xlsx"))
write_xlsx(DEG_down_KEGG_list,
           path = file.path(output_dir, "DEG_down_KEGG_results.xlsx"))

message("demo_bulk analysis completed. Results in: ", output_dir)
