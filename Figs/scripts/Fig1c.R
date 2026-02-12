## Fig 1c | DC bulk common genes heatmap
suppressPackageStartupMessages({
  library(pheatmap)
  library(readr)
  library(dplyr)
})

infile <- file.path("Figs", "data", "DCbulk_common_heatmap_exp.csv")
stopifnot(file.exists(infile))

outdir <- file.path("Figs", "outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, "Fig1c_DCbulk_common_heatmap.pdf")


merged_df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

# gene ids
rownames(merged_df) <- merged_df[[1]]
sample_cols <- 3:11
sample_names <- colnames(merged_df)[sample_cols]

# heatmap matrix: rows = samples, cols = genes
mat <- t(as.matrix(merged_df[, sample_cols]))
mode(mat) <- "numeric"
rownames(mat) <- sample_names
colnames(mat) <- rownames(merged_df)


# annotation_row: for samples (rows of mat)
treat_vec <- c(rep("WTDC", 3), rep("VecDC", 3), rep("iCDC", 3))

annotation_row <- data.frame(
  Treat = factor(treat_vec, levels = c("WTDC", "VecDC", "iCDC"))
)
rownames(annotation_row) <- rownames(mat)

annotation_col <- data.frame(
  GeneClass = factor(merged_df[[2]])
)
rownames(annotation_col) <- colnames(mat)


ann_colors <- list(
  Treat = c(iCDC = "#700303", WTDC = "grey30", VecDC = "grey60"),
  GeneClass = c(
    "Chemotaxis" = "#a9d6e8",
    "Immune Activation" = "#2799c8",
    "Pro-inflammatory" = "#70bca5",
    "Anti-inflammatory" = "#f6c06e",
    "pro-Treg Differentiation" = "#d9a4c4"
  )
)
hm_colors <- colorRampPalette(c("#4575B4", "#FEFEC0", "#D73027"))(250)


index_break <- c(6, 11, 20, 28)
index_break <- index_break[index_break > 0 & index_break < ncol(mat)]
if (length(index_break) == 0) index_break <- NULL


grDevices::cairo_pdf(outfile, width = 12, height = 3.2)  
pheatmap(
  mat,
  scale = "column",                 
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "white",
  color = hm_colors,

  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,

  gaps_col = index_break,          
  show_rownames = TRUE,
  show_colnames = FALSE,           
  fontsize_row = 9,
  main = "Gene Expression, DC bulk, Common"
)
dev.off()

message("Saved: ", outfile)
