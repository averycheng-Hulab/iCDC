## Fig 1d | T bulk common genes heatmap
suppressPackageStartupMessages({
  library(pheatmap)
})

infile <- file.path("Figs", "data", "DCbulk_common_heatmap_exp.csv")
stopifnot(file.exists(infile))

outdir <- file.path("Figs", "outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, "Fig1d_Tbulk_common_heatmap.pdf")

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

# annotations
annotation_row <- data.frame(
  Treat = factor(
    c(rep("TA", 3), rep("TA_VecDC", 3), rep("TA_iCDC", 3)),
    levels = c("TA", "TA_VecDC", "TA_iCDC")
  )
)
rownames(annotation_row) <- rownames(mat)

annotation_col <- data.frame(
  GeneClass = factor(merged_df[[2]])
)
rownames(annotation_col) <- colnames(mat)

ann_colors <- list(
  Treat = c(TA_iCDC = "#700303", TA = "grey60", TA_VecDC = "grey30"),
  GeneClass = c(
    "Chemotaxis" = "#a9d6e8",
    "Immune Activation" = "#2799c8",
    "Pro-inflammatory" = "#70bca5"
  )
)

hm_colors <- colorRampPalette(c("#4575B4", "#FEFEC0", "#D73027"))(250)

index_break <- c(4, 12)

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
  main = "Gene Expression, T bulk, Common"
)
dev.off()

message("Saved: ", outfile)
