############################################################
# bulk_1-2.R
# Bulk RNA-seq DEG enrichment (T cell / DC bulk)
#
# - Input: Excel DEG tables (FPKM-based, provided by platform)
#          columns must contain at least:
#          - gene_name
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
# NOTE:
#   Expression cutoff (exp_cutoff) is used only here because the
#   platform-derived DEG tables are FPKM-based and may contain
#   very low-expression genes. See README_bulk.md for details.
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# Source shared utilities
source(file.path("bulk", "utils_bulk.R"))

############################################################
# 1. User settings
############################################################

# Input directory containing DEG Excel files (FPKM-based)
input_dir  <- file.path("bulk", "input", "bulk_1-2")
output_dir <- file.path("bulk", "output", "bulk_1-2")
plot_dir   <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Replace these with your actual DEG filenames
deg_files <- c(
  file.path(input_dir, "DEG_Tcell_example.xlsx"),
  file.path(input_dir, "DEG_DC_example.xlsx")
)

# Names corresponding to each comparison
deg_names <- c("Tcell_Condition_vs_Control", "DC_Condition_vs_Control")

# Check consistency
stopifnot(length(deg_files) == length(deg_names))

# Parameters
exp_cutoff  <- 35   # minimal expression cutoff (applied to 2 expression columns)
logfc_cutoff <- 1   # absolute log2FC threshold
ntop        <- 3000 # maximum number of genes for enrichment

# Indices of expression columns used for exp_cutoff; adjust to match your tables
# e.g. column 8 and 9 in your original script
expr_col_idx <- c(8, 9)

############################################################
# 2. Load DEG tables
############################################################

DEG_list <- lapply(deg_files, function(f) {
  message("Reading: ", f)
  as.data.frame(read_xlsx(f))
})
names(DEG_list) <- deg_names

############################################################
# 3. Storage lists
############################################################

DEG_up_list       <- list()
DEG_down_list     <- list()
DEG_up_GO_list    <- list()
DEG_down_GO_list  <- list()
DEG_up_KEGG_list  <- list()
DEG_down_KEGG_list<- list()

############################################################
# 4. Main loop: filter, split, enrich
############################################################

for (i in seq_along(DEG_list)) {
  comp_name <- deg_names[i]
  df <- DEG_list[[i]]

  message("Processing comparison: ", comp_name)

  # Sanity checks
  required_cols <- c("gene_name", "log2FoldChange", "padj")
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
    dplyr::pull(gene_name) %>%
    unique()

  # Downregulated genes
  down_genes <- df %>%
    dplyr::filter(log2FoldChange < -logfc_cutoff) %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene_name) %>%
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

  # Example plots (top 10 terms/pathways), saved per comparison
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
}

############################################################
# 5. Optional: intersections between comparisons
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
# 6. Save results as Excel
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

message("bulk_1-2 analysis completed. Results in: ", output_dir)
