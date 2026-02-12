suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(eulerr)
})

# =========================
# I/O
# =========================
in_csv  <- file.path("Figs", "data", "cellstat_treg_expanded.csv")
out_dir <- file.path("Figs", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(in_csv, stringsAsFactors = FALSE)

# -------------------------
# Common settings
# -------------------------
cols_fill <- c("Non-Expanded" = "#2C7FB8", "Expanded" = "#F1A340")
y_breaks  <- (1:4) * 0.005

# =========================
# Part 1: Expanded vs Non-Expanded (%)
# =========================
make_stacked_percent <- function(df, out_pdf) {
  
  df_use <- df %>%
    filter(CellType0 == "T-Reg",
           orig.ident %in% c("MI_heart", "iCDC_heart"))
  
  plot_df <- df_use %>%
    count(orig.ident, exp_binary, name = "n") %>%
    group_by(orig.ident) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(
      exp_binary = factor(exp_binary, levels = c("Expanded","Non-Expanded")),
      orig.ident = factor(orig.ident, levels = c("MI_heart", "iCDC_heart"))
    )
  
  p <- ggplot(plot_df, aes(x = orig.ident, y = percent, fill = exp_binary)) +
    geom_col(width = 0.75, color = "black") +
    geom_text(
      aes(label = ifelse(percent >= 3, sprintf("%.1f", percent), "")),
      position = position_stack(vjust = 0.5),
      size = 4
    ) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = cols_fill) +
    labs(x = NULL, y = "Percentage (%)", fill = NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  ggsave(out_pdf, p, width = 3.2, height = 2.6)
}


# =========================
# Part 2/3: Expanded clonotype normalized count (per sample)
# =========================
make_expanded_clonotype_bar <- function(df, sample_id, out_pdf, title_text,
                                        y_breaks = (1:4) * 0.005,
                                        y_lim = c(0, 0.02)) {
  
  df_use <- df %>%
    filter(
      orig.ident == sample_id,
      CellType0 == "T-Reg",
      !is.na(clonotype),
      clonotype != ""
    )
  
  total_cells <- nrow(df_use)
  
  plot_df <- df_use %>%
    count(clonotype, name = "n_cells") %>%
    filter(n_cells > 1) %>%
    mutate(RelativeCellCount = n_cells / total_cells) %>%
    arrange(desc(RelativeCellCount))
  
  top_m <- nrow(plot_df)
  clonotype_order <- plot_df$clonotype
  
  p <- ggplot(
    plot_df,
    aes(x = factor(clonotype, levels = clonotype_order),
        y = RelativeCellCount)
  ) +
    geom_col(width = 0.8, color = "white", fill = "#F1A340") +
    scale_y_continuous(
      breaks = y_breaks,
      limits = y_lim,
      expand = expansion(mult = c(0, 0.1))
    ) +
    labs(
      x = "Clonotype ID",
      y = "Normalized Cell Count",
      title = paste0(title_text, " ( n = ", top_m, " )")
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_text(color = "black", face = "bold", size = 15, angle = 60, vjust = 0.5),
      axis.text.y  = element_text(color = "black", face = "bold", size = 15),
      axis.title.y = element_text(color = "black", face = "bold", size = 15),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 18)
    )
    ggsave(out_pdf,set_panel_size(p,width = unit(0.3*top_m,"in"),height = unit(6,"in")),width = 15,height = 8)
  
}


# =========================
# Part 4: Venn (iCDC_blood vs iCDC_heart)
# =========================
df <- df %>%
  filter(CellType0 == "T-Reg",
         !is.na(clonotype), clonotype != "")

icdc_b <- df %>% filter(orig.ident == "iCDC_blood") %>% pull(clonotype) %>% unique()
icdc_h <- df %>% filter(orig.ident == "iCDC_heart") %>% pull(clonotype) %>% unique()

dat <- c(
  "iCDC_blood" = length(icdc_b),
  "iCDC_heart" = length(icdc_h),
  "iCDC_blood&iCDC_heart" = length(intersect(icdc_b, icdc_h))
)

pdf(file.path(out_dir, "Fig3j_part4_Venn_iCDC_Treg.pdf"), height = 7, width = 7, onefile = TRUE)
plot(
  euler(dat),
  fills = list(fill = c("#d9a4c4", "#f6c06e"), alpha = 0.5),
  labels = list(col = "black", fontsize = 15),
  quantities = c(length(icdc_b),length(icdc_h),length(intersect(icdc_b,icdc_h))),
  main = "ClonoType Overlap"
)
dev.off()

# =========================
# Run all parts
# =========================
make_stacked_percent(df, file.path(out_dir, "Fig3j_part1_TregExpandedPercent.pdf"))

make_expanded_clonotype_bar(
  df, sample_id = "MI_heart",
  out_pdf = file.path(out_dir, "Fig3j_part2_MIheart_ExpandedClonotypes.pdf"),
  title_text = "Expanded Clonotypes, MI Heart",
  y_breaks = y_breaks, y_lim = c(0, 0.02)
)

make_expanded_clonotype_bar(
  df, sample_id = "iCDC_heart",
  out_pdf = file.path(out_dir, "Fig3j_part3_iCDCheart_ExpandedClonotypes.pdf"),
  title_text = "Expanded Clonotypes, iCDC Heart",
  y_breaks = y_breaks, y_lim = c(0, 0.02)
)
