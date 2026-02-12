## Fig 3a | Heart tissue GO/KEGG bar (Car vs IR, Down/Up)
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(gground)
  library(ggprism)
})

infile <- file.path("Figs", "data", "hearttissue_gokegg_bar.xlsx")
stopifnot(file.exists(infile))

outdir <- file.path("Figs", "outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pal <- c("#4575B4", "#c46878", "#BC8F8F")

plot_gokegg_bar <- function(sheet, title, outfile) {
  frame_plot_all <- data.frame(read_xlsx(infile, sheet = sheet))
  frame_plot_all <- frame_plot_all[!is.na(frame_plot_all[, 2]), ]

  frame_plot_all$Enrich <- factor(frame_plot_all$Enrich, levels = rev(c("GO", "KEGG")))
  frame_plot_all <- frame_plot_all %>% arrange(Enrich, desc(pvalue))
  frame_plot_all$Description <- factor(frame_plot_all$Description, levels = frame_plot_all$Description)
  frame_plot_all <- frame_plot_all %>% rowid_to_column("index")

  width <- 1
  xaxis_max <- max(-log10(frame_plot_all$pvalue)) + 1

  rect.data <- frame_plot_all %>%
    group_by(Enrich) %>%
    reframe(n = n()) %>%
    ungroup() %>%
    mutate(
      xmin = -3 * width,
      xmax = -2 * width,
      ymax = cumsum(n),
      ymin = lag(ymax, default = 0) + 0.6,
      ymax = ymax + 0.4
    )

  p <- frame_plot_all %>%
    ggplot(aes(-log10(pvalue), y = index, fill = Enrich)) +
    geom_round_col(aes(y = Description), width = 0.6, alpha = 0.8) +
    geom_text(aes(x = 0.2, label = Description), hjust = 0, size = 5, fontface = "bold") +
    geom_point(aes(x = -width, size = Count), shape = 21) +
    geom_text(aes(x = -width, label = Count), fontface = "bold") +
    scale_size_continuous(name = "Count", range = c(13, 13)) +
    geom_round_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Enrich),
      data = rect.data,
      radius = unit(4, "mm"),
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = Enrich),
      data = rect.data,
      angle = 90,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_segment(
      aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
      linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    labs(y = NULL) +
    scale_fill_manual(name = "Category", values = pal) +
    scale_x_continuous(breaks = seq(0, xaxis_max, 2), expand = expansion(c(0, 0))) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text()
    ) +
    ggtitle(title)

  ggsave(outfile, p, width = 12, height = 9.5)
  message("Saved: ", outfile)
}

plot_gokegg_bar(
  sheet = 1,
  title = "Car vs IR, Down",
  outfile = file.path(outdir, "Fig3a_HeartBulk_CarIR_GO_KEGG_Down.pdf")
)

plot_gokegg_bar(
  sheet = 2,
  title = "Car vs IR, Up",
  outfile = file.path(outdir, "Fig3a_HeartBulk_CarIR_GO_KEGG_Up.pdf")
)
