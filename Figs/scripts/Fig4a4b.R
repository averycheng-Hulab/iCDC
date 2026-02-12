## Fig 4a-b | Mac/Mono/DC UMAP + Mac composition
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

infile <- file.path("Figs", "data", "Macmonodc_plot_input.csv")
stopifnot(file.exists(infile))

outdir <- file.path("Figs", "outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

outfile_umap <- file.path(outdir, "Fig4a_Macmonodc_umap_by_Subtype.pdf")
outfile_comp <- file.path(outdir, "Fig4b_Mac_cellstat.pdf")

df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

# clean cluster labels once
df <- df %>%
  mutate(
    cluster_label = str_replace_all(cluster_label, "[\u2010\u2011\u2012\u2013\u2212]", "-"),
    cluster_label = str_trim(cluster_label)
  )

## ---------------------------
## Fig 4a | UMAP colored by subtype
## ---------------------------
clusters_umap <- c(
  "Lyve1+ Mac","Ccr2-MHCII-hi Mac","Ccr2+MHCII-hi Mac","Ccl8+ Mac","Spp1+ Mac",
  "IFN Mac","Fn1+ Mac","Saa3+ Mac","Prolif Mac","Ly6c2-hi Mono","Ly6c2-lo Mono",
  "cDC1","cDC2","pDC","Migratory DC","Other"
)

cols_umap <- c(
  "#a9d6e8","#2799c8","#2f6586","#70bca5","#82cf6c","#41af2e",
  "#f6c06e","#ff9400","#ff8f4f","#ff777f","#f3002b",
  "#d9a4c4","#9b71bb","#b59aa4","#f3dd84","#e2dfe0"
)

pal_umap <- setNames(cols_umap, clusters_umap)

df_umap <- df %>%
  mutate(cluster_label = factor(cluster_label, levels = clusters_umap))

p_umap <- ggplot(df_umap, aes(x = umap_1, y = umap_2, color = cluster_label)) +
  geom_point(size = 0.15, alpha = 0.8) +
  coord_equal() +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "right") +
  scale_color_manual(values = pal_umap, breaks = clusters_umap, drop = FALSE) +
  labs(x = "UMAP_1", y = "UMAP_2", color = NULL)

ggsave(outfile_umap, p_umap, width = 6.8, height = 5.2, useDingbats = FALSE)

## ---------------------------
## Fig 4b | Mac composition (stacked bar)
## ---------------------------
clusters_mac <- c(
  "Lyve1+ Mac","Ccr2-MHCII-hi Mac","Ccr2+MHCII-hi Mac","Ccl8+ Mac","Spp1+ Mac",
  "IFN Mac","Fn1+ Mac","Saa3+ Mac","Prolif Mac"
)

cols_mac <- c(
  "#a9d6e8","#2799c8","#2f6586","#70bca5","#82cf6c","#41af2e",
  "#f6c06e","#ff9400","#ff8f4f"
)

pal_mac <- setNames(cols_mac, clusters_mac)

df_mac <- df %>%
  filter(cluster_label0 == "Mac") %>%
  mutate(
    cluster_label = factor(cluster_label, levels = clusters_mac),
    Treatment = factor(Treatment, levels = c("Sham", "IR", "DC"))
  )

plot_df <- df_mac %>%
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
  scale_fill_manual(values = pal_mac, breaks = clusters_mac, drop = FALSE) +
  theme_bw() +
  labs(x = NULL, y = "Percentage (%)", fill = NULL, title = "Mac composition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(outfile_comp, p_comp, width = 7, height = 6, useDingbats = FALSE)

message("Saved: ", outfile_umap)
message("Saved: ", outfile_comp)
