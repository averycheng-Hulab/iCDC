############################################################
# bulk_3.R
# Bulk RNA-seq DEG enrichment (whole-heart bulk)
#
# - Input: CSV DEG tables generated from your DESeq2 pipeline
#          Expected columns include:
#          geneid_symbol, baseMean, log2FoldChange, lfcSE,
#          stat, pvalue, padj
#          plus TPM columns: S1...S6, V50, V52, V57, V65
#
# - Output: DEG up/down tables, GO/KEGG enrichment,
#           intersection & trend sets (as in your original script)
#
# - Dataset:
#   bulk_3: whole-heart bulk RNA-seq (TPM-based DESeq2)
#
# NOTE:
#   Here we do NOT use an expression (TPM) cutoff, because:
#   - DESeq2 starts from raw counts;
#   - The DEG tables are already filtered at the differential
#     expression step (p-value / FDR thresholds).
#   This rationale is described in README_bulk.md.
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(writexl)
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# Source shared utilities
source(file.path("bulk", "utils_bulk.R"))

############################################################
# 1. User settings
############################################################

input_dir  <- file.path("bulk", "input", "bulk_3")
output_dir <- file.path("bulk", "output", "bulk_3")
plot_dir   <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# DEG files (TPM-based, from your own DESeq2 pipeline)
deg_files <- c(
  "IR_vs_CAR.csv",   # originally 20250524_IR_Car_deseq2Results_2_tpm.csv
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

############################################################
# 2. Load DEG tables
############################################################

DEG_list <- lapply(deg_paths, function(f) {
  message("Reading: ", f)
  read.csv(f, header = TRUE, stringsAsFactors = FALSE)
})

names(DEG_list) <- contrast_names

############################################################
# 3. Storage lists
############################################################

DEG_up_list        <- list()
DEG_down_list      <- list()
DEG_up_GO_list     <- list()
DEG_down_GO_list   <- list()
DEG_up_KEGG_list   <- list()
DEG_down_KEGG_list <- list()

############################################################
# 4. Main loop: filter by log2FC & p-value, enrich
############################################################

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
}

############################################################
# 5. Intersections & trend gene sets (as in original script)
############################################################
# Here we follow your logic:
# - intersection of up DEGs between IR_vs_CAR and Vector_vs_CAR
# - intersection of up DEGs between Sham_vs_IR and Sham_vs_Vector
# - same for down DEGs
# - trend sets: up in (Car comparisons) & down in (Sham comparisons), etc.
############################################################

# Helper to get symbol vector from stored tables
get_symbols <- function(df_list_entry) {
  if (is.null(df_list_entry) || nrow(df_list_entry) == 0) {
    return(character(0))
  }
  gsub(" ", "", df_list_entry$geneid_symbol)
}

# Up DEGs intersection: IR_vs_CAR & Vector_vs_CAR
genes_up_12 <- intersect(
  get_symbols(DEG_up_list[["IR_vs_CAR"]]),
  get_symbols(DEG_up_list[["Vector_vs_CAR"]])
)

# Up DEGs intersection: Sham_vs_IR & Sham_vs_Vector
genes_up_34 <- intersect(
  get_symbols(DEG_up_list[["Sham_vs_IR"]]),
  get_symbols(DEG_up_list[["Sham_vs_Vector"]])
)

# Down DEGs intersection: IR_vs_CAR & Vector_vs_CAR
genes_down_12 <- intersect(
  get_symbols(DEG_down_list[["IR_vs_CAR"]]),
  get_symbols(DEG_down_list[["Vector_vs_CAR"]])
)

# Down DEGs intersection: Sham_vs_IR & Sham_vs_Vector
genes_down_34 <- intersect(
  get_symbols(DEG_down_list[["Sham_vs_IR"]]),
  get_symbols(DEG_down_list[["Sham_vs_Vector"]])
)

# Store intersection gene tables (for clarity)
DEG_up_list[["Up_IRvsCAR_and_VectorvsCAR"]]         <- data.frame(genes = genes_up_12)
DEG_up_list[["Up_ShamvsIR_and_ShamvsVector"]]       <- data.frame(genes = genes_up_34)
DEG_down_list[["Down_IRvsCAR_and_VectorvsCAR"]]     <- data.frame(genes = genes_down_12)
DEG_down_list[["Down_ShamvsIR_and_ShamvsVector"]]   <- data.frame(genes = genes_down_34)

# Enrich these intersection sets
DEG_up_GO_list[["Up_IRvsCAR_and_VectorvsCAR"]]       <- enrich_GO_KEGG_mouse(genes_up_12,   ntop)$GO
DEG_up_KEGG_list[["Up_IRvsCAR_and_VectorvsCAR"]]     <- enrich_GO_KEGG_mouse(genes_up_12,   ntop)$KEGG
DEG_up_GO_list[["Up_ShamvsIR_and_ShamvsVector"]]     <- enrich_GO_KEGG_mouse(genes_up_34,   ntop)$GO
DEG_up_KEGG_list[["Up_ShamvsIR_and_ShamvsVector"]]   <- enrich_GO_KEGG_mouse(genes_up_34,   ntop)$KEGG

DEG_down_GO_list[["Down_IRvsCAR_and_VectorvsCAR"]]   <- enrich_GO_KEGG_mouse(genes_down_12, ntop)$GO
DEG_down_KEGG_list[["Down_IRvsCAR_and_VectorvsCAR"]] <- enrich_GO_KEGG_mouse(genes_down_12, ntop)$KEGG
DEG_down_GO_list[["Down_ShamvsIR_and_ShamvsVector"]] <- enrich_GO_KEGG_mouse(genes_down_34, ntop)$GO
DEG_down_KEGG_list[["Down_ShamvsIR_and_ShamvsVector"]]<- enrich_GO_KEGG_mouse(genes_down_34, ntop)$KEGG

############################################################
# 6. Trend gene sets (Up→Down and Down→Up)
############################################################
# Example interpretation (you can phrase biologically in manuscript):
# - Trend_UpInCAR_DownInSham: genes consistently up in CAR-related
#   contrasts and down in Sham-related contrasts.
############################################################

DEG_trend_list       <- list()
DEG_trend_GO_KEGG    <- list()

# Up in CAR-related (intersection 1&2), Down in Sham-related (intersection 3&4)
genes_trend_updown <- intersect(genes_up_12, genes_down_34)
DEG_trend_list[["UpIn_CAR_downIn_Sham"]] <- data.frame(genes = genes_trend_updown)

enrich_updown <- enrich_GO_KEGG_mouse(genes_trend_updown, ntop = ntop)
DEG_trend_GO_KEGG[["UpIn_CAR_downIn_Sham_GO"]]   <- enrich_updown$GO
DEG_trend_GO_KEGG[["UpIn_CAR_downIn_Sham_KEGG"]] <- enrich_updown$KEGG

# Down in CAR-related, Up in Sham-related
genes_trend_downup <- intersect(genes_down_12, genes_up_34)
DEG_trend_list[["DownIn_CAR_upIn_Sham"]] <- data.frame(genes = genes_trend_downup)

enrich_downup <- enrich_GO_KEGG_mouse(genes_trend_downup, ntop = ntop)
DEG_trend_GO_KEGG[["DownIn_CAR_upIn_Sham_GO"]]   <- enrich_downup$GO
DEG_trend_GO_KEGG[["DownIn_CAR_upIn_Sham_KEGG"]] <- enrich_downup$KEGG

############################################################
# 7. Save results
############################################################

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

# Trend sets
write_xlsx(DEG_trend_list,
           path = file.path(output_dir, "DEG_trend_gene_sets.xlsx"))
write_xlsx(DEG_trend_GO_KEGG,
           path = file.path(output_dir, "DEG_trend_GO_KEGG_results.xlsx"))

# Save R objects if needed for further downstream plotting
save(DEG_trend_GO_KEGG,
     file = file.path(output_dir, "DEG_trend_GO_KEGG_list.RData"))
save(DEG_trend_list,
     file = file.path(output_dir, "DEG_trend_list.RData"))
save(DEG_up_list,
     file = file.path(output_dir, "DEG_up_list.RData"))
save(DEG_down_list,
     file = file.path(output_dir, "DEG_down_list.RData"))
save(DEG_up_GO_list,
     file = file.path(output_dir, "DEG_up_GO_list.RData"))
save(DEG_down_GO_list,
     file = file.path(output_dir, "DEG_down_GO_list.RData"))
save(DEG_up_KEGG_list,
     file = file.path(output_dir, "DEG_up_KEGG_list.RData"))
save(DEG_down_KEGG_list,
     file = file.path(output_dir, "DEG_down_KEGG_list.RData"))

message("bulk_3 analysis completed. Results in: ", output_dir)
