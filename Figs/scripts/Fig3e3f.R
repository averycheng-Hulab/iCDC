## Fig 3e-f | TNK UMAP + T subset composition
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

infile <- file.path("Figs", "data", "Tnk_plot_input.csv")
stopifnot(file.exists(infile))

outdir <- file.path("Figs", "outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

outfile_umap <- file.path(outdir, "Fig3e_TNK_umap_by_Subtype.pdf")
outfile_comp <- file.path(outdir, "Fig3f_T_cellstat.pdf")

df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

# clean cluster labels once (keep consistent across plots)
df <- df %>%
  mutate(
    cluster_label = str_replace_all(cluster_label, "[\u2010\u2011\u2012\u2013\u2212]", "-"),
    cluster_label = str_trim(cluster_label)
  )

## ---------------------------
## 1) UMAP colored by subtype
## ---------------------------
clusters_umap <- c(
  "T-Helper","Cytotoxic","gdT","Treg","NaiveT","TCM","T-DN","NKT","NK","ILC","Other"
)
cols_umap <- c(
  "#A1C5D8","#246DA0","#A8CC85","#3A9338","#EEB66E","#E47920","#E89191","#C9241F","#C2ADCD","#603888","#F1EE9B"
)
pal_umap <- setNames(cols_umap, clusters_umap)

df <- df %>%
  mutate(cluster_label = factor(cluster_label, levels = clusters_umap))

p_umap <- ggplot(df, aes(x = umap_1, y = umap_2, color = cluster_label)) +
  geom_point(size = 0.15, alpha = 0.8) +
  coord_equal() +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "right"
  ) +
  scale_color_manual(values = pal_umap, breaks = clusters_umap, drop = FALSE) +
  labs(x = "UMAP_1", y = "UMAP_2", color = NULL)

ggsave(outfile_umap, p_umap, width = 6.8, height = 5.2)

## ---------------------------
## 2) T subset composition (stacked bar)
## ---------------------------
clusters_t <- c("T-Helper","Cytotoxic","gdT","Treg","NaiveT","TCM","T-DN","NKT")
cols_t <- c("#A1C5D8","#246DA0","#A8CC85","#3A9338","#EEB66E","#E47920","#E89191","#C9241F")
pal_t <- setNames(cols_t, clusters_t)

df_tnk <- df %>%
  filter(cluster_label0 == "T") %>%
  mutate(
    cluster_label = factor(as.character(cluster_label), levels = clusters_t),
    Treatment = factor(Treatment, levels = c("Sham", "IR", "DC"))
  )

plot_df <- df_tnk %>%
  count(Treatment, cluster_label, name = "N") %>%
  group_by(Treatment) %>%
  mutate(percent = N / sum(N) * 100) %>%
  ungroup()

p_comp <- ggplot(plot_df, aes(x = Treatment, y = percent, fill = cluster_label)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = ifelse(percent >= 6, sprintf("%.1f", percent), "")),
    position = position_stack(vjust = 0.5),
    size = 5
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = pal_t, breaks = clusters_t, drop = FALSE) +
  theme_bw() +
  labs(x = NULL, y = "Percentage (%)", fill = NULL, title = "T composition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(outfile_comp, p_comp, width = 7, height = 6)

message("Saved: ", outfile_umap)
message("Saved: ", outfile_comp)
