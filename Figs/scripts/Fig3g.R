## Fig 3g | scRNA T common genes heatmap
suppressPackageStartupMessages({
  library(pheatmap)
})

infile <- file.path("Figs", "data", "sc_T_common_heatmap_exp.csv")
stopifnot(file.exists(infile))

outdir <- file.path("Figs", "outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, "Fig3g_scT_common_heatmap.pdf")

merged_df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

# gene ids
rownames(merged_df) <- merged_df[[1]]
sample_cols <- 3:13
sample_names <- colnames(merged_df)[sample_cols]

# heatmap matrix: rows = samples, cols = genes
mat <- t(as.matrix(merged_df[, sample_cols]))
mode(mat) <- "numeric"
rownames(mat) <- sample_names
colnames(mat) <- rownames(merged_df)

annotation_row <- data.frame(
  Treat = factor(
    c(rep("Sham", 3), rep("IR", 4), rep("iCDC", 4)),
    levels = c("Sham", "IR", "iCDC")
  )
)
rownames(annotation_row) <- rownames(mat)

annotation_col <- data.frame(
  GeneClass = factor(merged_df[[2]])
)
rownames(annotation_col) <- colnames(mat)

ann_colors <- list(
  Treat = c(iCDC = "#700303", Sham = "grey30", IR = "grey60"),
  GeneClass = c(
    "Cell Migration" = "#a9d6e8",
    "Pro-inflammation" = "#2799c8",
    "Th1/2/17/T-Toxic  Cell Differentiation" = "#d9a4c4"
  )
)

hm_colors <- colorRampPalette(c("#4575B4", "#FEFEC0", "#D73027"))(250)

index_break <- c(8, 12)

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
  main = "Gene Expression, SC, T, Common"
)
dev.off()

message("Saved: ", outfile)
