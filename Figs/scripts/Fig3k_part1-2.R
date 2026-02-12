suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# =========================
# I/O
# =========================
in_dir  <- file.path("Figs", "data")
out_dir <- file.path("Figs", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# Helper
# =========================
make_violin_panel <- function(in_csv, out_pdf, gene_list,
                              group_levels = c("MI_heart", "iCDC_heart")) {

  df <- read.csv(in_csv, stringsAsFactors = FALSE) %>%
    mutate(
      orig.ident = factor(orig.ident, levels = group_levels),
      Gene = factor(Gene, levels = gene_list)
    )

  p <- ggplot(df, aes(x = Gene, y = Expression, fill = orig.ident)) +
    geom_violin(scale = "width", position = position_dodge(0.8),
                trim = TRUE, alpha = 0.7, width = 0.7) +
    geom_boxplot(width = 0.15, position = position_dodge(0.8),
                 outlier.shape = NA, outlier.alpha = 0.3) +
    stat_summary(fun = "median", geom = "point",
                 position = position_dodge(0.8),
                 size = 1, color = "white") +
    scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
    labs(x = NULL, y = "Expression Level", fill = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
      legend.position = "top"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

  ggsave(out_pdf, p, width = 8, height = 3, dpi = 300)
}

# =========================
# Fig3k
# =========================
genes_treg <- c(
  "Ccr8","Maf","Tnfrsf4","Tnfrsf9","Tgfb1","Tnfrsf18",
  "Ets1","Gimap3","Gimap5","Gpx1","Gpx4","Il2rg","Lgals1","Il2rb",
  "Nfkb1","Cox5a","Txn1","Prdx2","Prdx1"
)

genes_treg_exp <- c(
  "Sipa1","Trbc2","Trbv29","Gzmb","Lime1","Selenok",
  "Il31ra","Aven","Cct6a","Ppm1h","Wrn"
)

make_violin_panel(
  in_csv  = file.path(in_dir, "violin_treg.csv"),
  out_pdf = file.path(out_dir, "Fig3k_part1_Treg_genes_violin.pdf"),
  gene_list = genes_treg
)

make_violin_panel(
  in_csv  = file.path(in_dir, "violin_treg_expanded.csv"),
  out_pdf = file.path(out_dir, "Fig3k_part2_ExpandedTreg_genes_violin.pdf"),
  gene_list = genes_treg_exp
)
