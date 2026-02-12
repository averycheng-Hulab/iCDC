## Fig 4j-k-l | Fib UMAP + Fib composition + Fib ECM Regulator Score
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggpubr)
})

infile <- file.path("Figs", "data", "Fib_plot_input.csv")
stopifnot(file.exists(infile))

outdir <- file.path("Figs", "outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

outfile_umap <- file.path(outdir, "Fig4j_Fib_umap_by_Subtype.pdf")
outfile_comp <- file.path(outdir, "Fig4k_Fib_cellstat.pdf")
outfile_ecm  <- file.path(outdir, "Fig4l_Fib_ECMRegScore_violin.pdf")

df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

# clean labels once
df <- df %>%
  mutate(
    cluster_label = str_replace_all(cluster_label, "[\u2010\u2011\u2012\u2013\u2212]", "-"),
    cluster_label = str_trim(cluster_label)
  )

clusters <- c("F-SL","F-SH","F-Myo","F-Act","F-IFNs","F-IR")
cols <- c("#A1C5D8","#246DA0","#A8CC85","#3A9338","#EEB66E","#E47920")
pal_fib <- setNames(cols, clusters)

## ---------------------------
## Fig 4j | UMAP colored by subtype
## ---------------------------
df_umap <- df %>%
  mutate(cluster_label = factor(cluster_label, levels = clusters))

p_umap <- ggplot(df_umap, aes(x = umap_1, y = umap_2, color = cluster_label)) +
  geom_point(size = 0.15, alpha = 0.8) +
  coord_equal() +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "right") +
  scale_color_manual(values = pal_fib, breaks = clusters, drop = FALSE) +
  labs(x = "UMAP_1", y = "UMAP_2", color = NULL)

ggsave(outfile_umap, p_umap, width = 6.8, height = 5.2, useDingbats = FALSE)

## ---------------------------
## Fig 4k | Fib composition (stacked bar)
## ---------------------------
df_fib <- df %>%
  mutate(
    cluster_label = factor(cluster_label, levels = clusters),
    Treatment = factor(Treatment, levels = c("Sham", "IR", "Vector", "DC"))
  )

plot_df <- df_fib %>%
  count(Treatment, cluster_label, name = "n") %>%
  group_by(Treatment) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

p_comp <- ggplot(plot_df, aes(x = Treatment, y = percent, fill = cluster_label)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = ifelse(percent >= 6, sprintf("%.1f", percent), "")),
    position = position_stack(vjust = 0.5),
    size = 5
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = pal_fib, breaks = clusters, drop = FALSE) +
  theme_bw() +
  labs(x = NULL, y = "Percentage (%)", fill = NULL, title = "Fib composition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(outfile_comp, p_comp, width = 7, height = 6, useDingbats = FALSE)

## ---------------------------
## Fig 4l | ECM Regulator Score (violin)
## ---------------------------
score_col <- if ("ECM_Regulator_Score" %in% colnames(df)) {
  "ECM_Regulator_Score"
} else if ("ECM Regulator Score" %in% colnames(df)) {
  "ECM Regulator Score"
} else {
  stop("找不到 ECM score 列：请确认 input.csv 里有 ECM_Regulator_Score 或 `ECM Regulator Score`")
}

plot_df2 <- df %>%
  dplyr::select(Treatment, score = all_of(score_col)) %>%
  dplyr::filter(!is.na(Treatment), !is.na(score))

plot_df2$Treatment <- factor(plot_df2$Treatment, levels = c("Sham", "IR", "Vector", "DC"))

comparisons <- list(
  c("Sham", "IR"),
  c("Sham", "Vector"),
  c("IR", "DC"),
  c("Vector", "DC")
)

cols4 <- c(
  "Sham"   = "#A1C5D8",
  "IR"     = "#246DA0",
  "Vector" = "#A8CC85",
  "DC"     = "#3A9338"
)

p_ecm <- ggplot(plot_df2, aes(x = Treatment, y = score, fill = Treatment)) +
  geom_violin(trim = TRUE, color = "grey30", linewidth = 0.35) +
  geom_boxplot(width = 0.18, outlier.shape = NA, color = "black", fill = "black", alpha = 0.18) +
  stat_summary(fun = median, geom = "point", size = 2.2, color = "white") +
  ggpubr::stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    hide.ns = TRUE,
    tip.length = 0.01,
    bracket.size = 0.25,
    step.increase = 0.10,
    size = 4
  ) +
  scale_fill_manual(values = cols4, breaks = names(cols4), drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Expression level", title = "ECM Regulator Score, Fibroblast") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(outfile_ecm, p_ecm, width = 3.2, height = 3.2, useDingbats = FALSE)

message("Saved: ", outfile_umap)
message("Saved: ", outfile_comp)
message("Saved: ", outfile_ecm)
