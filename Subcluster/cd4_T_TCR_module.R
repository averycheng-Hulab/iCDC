###############################################################################
# cd4_T_TCR_module.R
#
# - Input : r_objects/seuobj_cd4_final.Rds
# - Output: r_objects/seuobj_cd4_tcr_final.Rds
#
# - TCR files (standardized file names):
#   filtered_contig_annotations_MI_Blood.csv
#   filtered_contig_annotations_MI_Heart.csv
#   filtered_contig_annotations_icDC_Blood.csv
#   filtered_contig_annotations_icDC_Heart.csv
#
# - Steps:
#   (1) Load CD4 Seurat object
#   (2) Combine TCR contigs
#   (3) Map TCR info into Seurat (scRepertoire)
#   (4) Subset Treg
#   (5) Expanded vs Non-expanded definition
#   (6) Expansion barplots
#   (7) Save Seurat object
###############################################################################

## ---------------------------------------------------------------------------
## 0. Packages and load utilities
## ---------------------------------------------------------------------------

required_pkgs <- c(
  "Seurat","dplyr","ggplot2","patchwork","stringr",
  "scRepertoire"
)
missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse=", "))
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(scRepertoire)

source("utils_scRNA.R")   # ★ 使用统一的 plot_theme_common / save_plot

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)

dir.create("output_plots", showWarnings = FALSE)
dir.create("r_objects", showWarnings = FALSE)


## ---------------------------------------------------------------------------
## 1. Load CD4 final Seurat object
## ---------------------------------------------------------------------------

seuobj_cd4_final <- readRDS("r_objects/seuobj_cd4_final.Rds")


## ---------------------------------------------------------------------------
## 2. Load contig files
## ---------------------------------------------------------------------------

contigs_list <- list(
  read.csv("filtered_contig_annotations_MI_Blood.csv",  stringsAsFactors = FALSE),
  read.csv("filtered_contig_annotations_MI_Heart.csv",  stringsAsFactors = FALSE),
  read.csv("filtered_contig_annotations_icDC_Blood.csv", stringsAsFactors = FALSE),
  read.csv("filtered_contig_annotations_icDC_Heart.csv", stringsAsFactors = FALSE)
)

sample_names <- c("MI_Blood","MI_Heart","icDC_Blood","icDC_Heart")


## ---------------------------------------------------------------------------
## 3. Combine TCR with Seurat
## ---------------------------------------------------------------------------

combined_TCR <- combineTCR(
  contigs_list,
  samples = sample_names,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)

sce <- combineExpression(
  combined_TCR,
  seuobj_cd4_final,
  cloneCall = "aa",
  proportion = TRUE
)


## ---------------------------------------------------------------------------
## 4. Clonotype ID assignment
## ---------------------------------------------------------------------------

CTaa <- unique(na.omit(sce$CTaa))
clonotype_factor <- factor(paste0("CL", seq_along(CTaa)))
clonemap <- data.frame(CTaa = CTaa, clonotype = clonotype_factor)

meta <- sce@meta.data
meta2 <- dplyr::left_join(meta, clonemap, by = "CTaa")
rownames(meta2) <- rownames(meta)
sce@meta.data <- meta2


## ---------------------------------------------------------------------------
## 5. Subset to T-reg cells
## ---------------------------------------------------------------------------

sce_treg <- subset(sce, subset = CellType1 == "T-Reg")

sce_MI_Blood   <- subset(sce_treg, subset = orig.ident == "MI_Blood")
sce_MI_Heart   <- subset(sce_treg, subset = orig.ident == "MI_Heart")
sce_icDC_Blood <- subset(sce_treg, subset = orig.ident == "icDC_Blood")
sce_icDC_Heart <- subset(sce_treg, subset = orig.ident == "icDC_Heart")


## ---------------------------------------------------------------------------
## 6. Function: annotate clonotype expansion
## ---------------------------------------------------------------------------

annotate_expansion <- function(seu) {
  meta <- seu@meta.data
  meta_non_na <- meta[!is.na(meta$clonotype), ]
  meta_non_na$barcode <- rownames(meta_non_na)
  meta_non_na$count <- 1

  clonotype_freq <- seu@meta.data %>%
    filter(!is.na(clonotype)) %>%
    group_by(clonotype) %>%
    summarise(barcode_count = n()) %>%
    arrange(desc(barcode_count))

  merged <- merge(meta_non_na, clonotype_freq, by="clonotype", all.x=TRUE)
  merged <- merged %>% arrange(desc(barcode_count))
  merged$RelativeCellCount <- merged$count / length(merged$barcode_count)
  merged$exp_binary <- ifelse(merged$barcode_count > 1, "Expanded", "Non-Expanded")

  rownames(merged) <- merged$barcode
  merged <- merged[rownames(seu@meta.data), ]

  seu@meta.data <- merged
  return(seu)
}

sce_MI_Blood   <- annotate_expansion(sce_MI_Blood)
sce_MI_Heart   <- annotate_expansion(sce_MI_Heart)
sce_icDC_Blood <- annotate_expansion(sce_icDC_Blood)
sce_icDC_Heart <- annotate_expansion(sce_icDC_Heart)


## ---------------------------------------------------------------------------
## 7. Merge MI_Heart + icDC_Heart (Heart panel)
## ---------------------------------------------------------------------------

sce_heart <- merge(sce_MI_Heart, sce_icDC_Heart)


## ---------------------------------------------------------------------------
## 8. Expansion barplot (Heart)
## ---------------------------------------------------------------------------

p_expansion <- ggplot(
  sce_heart@meta.data,
  aes(x = orig.ident, fill = exp_binary)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title="Expanded Clonotypes (Heart)", x="", y="Fraction") +
  plot_theme_common()

save_plot(
  p_expansion,
  "cd4T_TCR_expansion_barplot.pdf",
  width = 6, height = 5
)


## ---------------------------------------------------------------------------
## 9. Expanded clonotype ranking (MI Heart)
## ---------------------------------------------------------------------------

meta_heart <- sce_MI_Heart@meta.data
meta_exp   <- meta_heart[!is.na(meta_heart$clonotype), ]
meta_exp$barcode <- rownames(meta_exp)
meta_exp$count <- 1

clonotype_freq <- sce_MI_Heart@meta.data %>%
  filter(!is.na(clonotype)) %>%
  group_by(clonotype) %>%
  summarise(barcode_count = n()) %>%
  arrange(desc(barcode_count))

merged_df <- merge(meta_exp, clonotype_freq, by="clonotype", all.x=TRUE)
merged_df <- merged_df %>% arrange(desc(barcode_count))
merged_df$RelativeCellCount <- merged_df$count / length(merged_df$barcode_count)

top_n <- sum(merged_df$barcode_count > 1)
top_m <- length(unique(merged_df[1:top_n, ]$clonotype))

frame_plot <- head(merged_df, top_n)
frame_plot$clonotype <- factor(frame_plot$clonotype, levels = unique(frame_plot$clonotype))

cols <- c(
  "T-CM/Naive"="#a9d6e8","T-Reg"="#2799c8",
  "T-EarlyAct"="#70bca5","T-Isg"="#82cf6c",
  "T-Trbv1"="#41af2e","T-EM"="#f6c06e",
  "T-Cycling"="#ff9400"
)

p_clone_rank <- ggplot(
  frame_plot,
  aes(x = clonotype, y = RelativeCellCount, fill = CellType1)
) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = cols) +
  labs(
    x = "Clonotype ID",
    y = "Cell Ratio of Expanded Clonotype (Treg)",
    title = paste("Expanded Clonotypes (n =", top_m, "), MI Heart")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title  = element_text(face = "bold"),
    plot.title  = element_text(face = "bold", hjust = 0.5)
  )

save_plot(
  p_clone_rank,
  "cd4T_TCR_cloneRanking_MIHeart.pdf",
  width = 12, height = 6
)


## ---------------------------------------------------------------------------
## 10. Save final Seurat object
## ---------------------------------------------------------------------------

seuobj_cd4_tcr_final <- sce
saveRDS(seuobj_cd4_tcr_final, "r_objects/seuobj_cd4_tcr_final.Rds")

message("cd4_T_TCR_module finished → r_objects/seuobj_cd4_tcr_final.Rds")

###############################################################################
# End of cd4_T_TCR_module.R
###############################################################################
