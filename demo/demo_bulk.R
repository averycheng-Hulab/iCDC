if (!dir.exists("demo") || !dir.exists("bulk")) stop("Please run from repository root.")

suppressPackageStartupMessages({
  library(readr)
  library(writexl)
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
})

# Source shared utilities
source(file.path("bulk","utils","utils_bulk.R"))

input_dir  <- file.path("demo","bulk_input")
output_dir <- file.path("demo","bulk_output")
plot_dir   <- output_dir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# DEG files 
deg_files <- c(
  "IR_vs_CAR.csv",  
  "Vector_vs_CAR.csv",
  "Sham_vs_IR.csv",
  "Sham_vs_Vector.csv"
)


contrast_names <- c("IR_vs_CAR", "Vector_vs_CAR", "Sham_vs_IR", "Sham_vs_Vector")

deg_paths <- file.path(input_dir, deg_files)
names(deg_paths) <- contrast_names

# Filtering parameters for DEGs
logfc_cutoff  <- 1      # |log2FC| > 1
pvalue_cutoff <- 0.01   # p-value threshold
ntop          <- 3000


## 1. Load DEG tables

DEG_list <- lapply(deg_paths, function(f) {
  message("Reading: ", f)
  read.csv(f, header = TRUE, stringsAsFactors = FALSE)
})
names(DEG_list) <- contrast_names


## 3. Storage lists

DEG_up_list        <- list()
DEG_down_list      <- list()
DEG_up_GO_list     <- list()
DEG_down_GO_list   <- list()
DEG_up_KEGG_list   <- list()
DEG_down_KEGG_list <- list()


## 2. Main loop: filter by log2FC & p-value, enrich

for (comp_name in names(DEG_list)) {
  df <- DEG_list[[comp_name]]
  message("Processing comparison: ", comp_name)

  # Sanity checks
  required_cols <- c("geneid_symbol", "log2FoldChange", "pvalue", "padj")
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
  up_genes <- gsub(" ", "", up_df$geneid_symbol)
  down_genes <- gsub(" ", "", down_df$geneid_symbol)

  # Enrichment
  up_enrich   <- enrich_GO_KEGG_mouse(up_genes,   ntop = ntop)
  down_enrich <- enrich_GO_KEGG_mouse(down_genes, ntop = ntop)

  DEG_up_GO_list[[comp_name]]     <- up_enrich$GO
  DEG_up_KEGG_list[[comp_name]]   <- up_enrich$KEGG
  DEG_down_GO_list[[comp_name]]   <- down_enrich$GO
  DEG_down_KEGG_list[[comp_name]] <- down_enrich$KEGG

  # Example plots (upregulated GO/KEGG)
  p_go   <- plot_GO_bar(up_enrich$GO,   n = 10,
                        title = paste("Top GO BP (Up) -", comp_name))
  p_kegg <- plot_KEGG_dot(up_enrich$KEGG, n = 10,
                          title = paste("Top KEGG (Up) -", comp_name))

  if (!is.null(p_go)) {
    ggsave(file.path(plot_dir, paste0("GO_up_", comp_name, ".pdf")),
           p_go, width = 7, height = 5)
  }
  if (!is.null(p_kegg)) {
    ggsave(file.path(plot_dir, paste0("KEGG_up_", comp_name, ".pdf")),
           p_kegg, width = 7, height = 5)
  }

  ###
  p_go2   <- plot_GO_bar(down_enrich$GO,   n = 10,
                        title = paste("Top GO BP (Down) -", comp_name))
  p_kegg2 <- plot_KEGG_dot(down_enrich$KEGG, n = 10,
                          title = paste("Top KEGG (Down) -", comp_name))

  if (!is.null(p_go2)) {
    ggsave(file.path(plot_dir, paste0("GO_down_", comp_name, ".pdf")),
           p_go2, width = 7, height = 5)
  }
  if (!is.null(p_kegg2)) {
    ggsave(file.path(plot_dir, paste0("KEGG_down_", comp_name, ".pdf")),
           p_kegg2, width = 7, height = 5)
  }
}


## 3. Save results
# DEG tables (up/down)
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
