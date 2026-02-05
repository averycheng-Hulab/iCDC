###############################################################################
# enrich_scRNA_cd4.R
#
# - Input : r_objects/seuobj_cd4_final.Rds
# - Contrast: MI vs iCDC (Treatment)
# - Target: one CD4 T-cell subset (e.g. T-Reg)
#
# - Output files (using cd4_subset name, e.g., "T-Reg"):
#   GoKegg_cd4/
#      T-Reg_feature_up.xlsx
#      T-Reg_feature_down.xlsx
#      T-Reg_GO_up.xlsx
#      T-Reg_GO_down.xlsx
#      T-Reg_KEGG_up.xlsx
#      T-Reg_KEGG_down.xlsx
#      T-Reg_GO_Up_bubble.pdf
#      T-Reg_GO_Down_bubble.pdf
#      T-Reg_KEGG_Up_bubble.pdf
#      T-Reg_KEGG_Down_bubble.pdf
###############################################################################

## ---------------------------------------------------------------------------
## 0. Packages, utils and options
## ---------------------------------------------------------------------------

required_pkgs <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "clusterProfiler",
  "DOSE",
  "org.Mm.eg.db",
  "writexl"
)

missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nPlease install them before running."
  )
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(writexl)

source("utils_scRNA.R")

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3)

dir.create("GoKegg_cd4", showWarnings = FALSE)


## ---------------------------------------------------------------------------
## 1. Parameters
## ---------------------------------------------------------------------------

# CD4 T-cell subset to analyze: "T-CM/Naive", "T-Reg", "T-EM", etc.
cd4_subset <- "T-Reg"

logfc_threshold <- 0.3
min_pct         <- 0.5
pval_cutoff     <- 0.01
ntop_go_genes   <- 100

ident_1      <- "MI"
ident_2      <- "iCDC"
group_by_var <- "Treatment"


## ---------------------------------------------------------------------------
## 2. Load Seurat object and subset the target CD4 subset
## ---------------------------------------------------------------------------

seuobj_cd4_final <- readRDS("r_objects/seuobj_cd4_final.Rds")

if (!group_by_var %in% colnames(seuobj_cd4_final@meta.data)) {
  stop("Column ", group_by_var, " not found in meta.data.")
}

if (!"CellType1" %in% colnames(seuobj_cd4_final@meta.data)) {
  stop("Column CellType1 not found in meta.data.")
}

seu_sub <- subset(seuobj_cd4_final, subset = CellType1 == cd4_subset)

if (ncol(seu_sub) == 0) {
  stop("No cells found for subset ", cd4_subset, ".")
}

if (!ident_1 %in% seu_sub@meta.data[[group_by_var]] ||
    !ident_2 %in% seu_sub@meta.data[[group_by_var]]) {
  stop(
    "Both ", ident_1, " and ", ident_2,
    " must be present in ", group_by_var, " for subset ", cd4_subset, "."
  )
}

if ("SCT" %in% Assays(seu_sub)) {
  DefaultAssay(seu_sub) <- "SCT"
}


## ---------------------------------------------------------------------------
## 3. Differential expression: MI vs iCDC
## ---------------------------------------------------------------------------

deg <- FindMarkers(
  object          = seu_sub,
  ident.1         = ident_1,
  ident.2         = ident_2,
  group.by        = group_by_var,
  assay           = "SCT",
  logfc.threshold = logfc_threshold,
  min.pct         = min_pct
)

deg$gene <- rownames(deg)

feature_up <- deg %>%
  filter(avg_log2FC > 0, p_val < pval_cutoff) %>%
  pull(gene) %>%
  unique()

feature_down <- deg %>%
  filter(avg_log2FC < 0, p_val < pval_cutoff) %>%
  pull(gene) %>%
  unique()

write_xlsx(
  data.frame(gene = feature_up),
  file.path("GoKegg_cd4", paste0(cd4_subset, "_feature_up.xlsx"))
)

write_xlsx(
  data.frame(gene = feature_down),
  file.path("GoKegg_cd4", paste0(cd4_subset, "_feature_down.xlsx"))
)


## ---------------------------------------------------------------------------
## 4. GO + KEGG enrichment helper
## ---------------------------------------------------------------------------

run_go_kegg <- function(gene_symbols) {
  genes_use <- gene_symbols[1:min(length(gene_symbols), ntop_go_genes)]
  genes_use <- genes_use[!is.na(genes_use)]

  if (length(genes_use) == 0) {
    return(list(GO = NULL, KEGG = NULL))
  }

  # GO enrichment
  ego <- enrichGO(
    gene          = genes_use,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    GO_df <- simplify(
      ego,
      cutoff     = 0.7,
      by         = "p.adjust",
      select_fun = min
    ) |>
      as.data.frame()
  } else {
    GO_df <- NULL
  }

  # KEGG enrichment: convert symbols to Entrez IDs
  genes_id <- bitr(
    genes_use,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  )

  if (!is.null(genes_id) && nrow(genes_id) > 0) {
    ekegg <- enrichKEGG(
      gene          = genes_id$ENTREZID,
      keyType       = "kegg",
      organism      = "mmu",
      pvalueCutoff  = 1,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1
    )

    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      kegg_df <- as.data.frame(ekegg@result)

      symbol_vec <- character(nrow(kegg_df))
      for (j in seq_len(nrow(kegg_df))) {
        gene_ids <- unlist(strsplit(kegg_df$geneID[j], "/"))
        gene2 <- bitr(
          gene_ids,
          fromType = "ENTREZID",
          toType   = "SYMBOL",
          OrgDb    = org.Mm.eg.db
        )
        symbol_vec[j] <- paste(gene2$SYMBOL, collapse = "/")
      }

      kegg_df$symbol_vec <- symbol_vec
      kegg_df$Description <- gsub(
        " - Mus musculus \\(house mouse\\)",
        "",
        kegg_df$Description
      )

      if ("category" %in% colnames(kegg_df)) {
        kegg_df <- kegg_df[kegg_df$category != "Human Diseases", ]
      }
    } else {
      kegg_df <- NULL
    }
  } else {
    kegg_df <- NULL
  }

  list(GO = GO_df, KEGG = kegg_df)
}


## ---------------------------------------------------------------------------
## 5. Run GO/KEGG enrichment for up and down genes
## ---------------------------------------------------------------------------

res_up    <- run_go_kegg(feature_up)
GO_up_df  <- res_up$GO
KEGG_up_df <- res_up$KEGG

res_down      <- run_go_kegg(feature_down)
GO_down_df    <- res_down$GO
KEGG_down_df  <- res_down$KEGG

if (!is.null(GO_up_df)) {
  write_xlsx(
    GO_up_df,
    file.path("GoKegg_cd4", paste0(cd4_subset, "_GO_up.xlsx"))
  )
}

if (!is.null(GO_down_df)) {
  write_xlsx(
    GO_down_df,
    file.path("GoKegg_cd4", paste0(cd4_subset, "_GO_down.xlsx"))
  )
}

if (!is.null(KEGG_up_df)) {
  write_xlsx(
    KEGG_up_df,
    file.path("GoKegg_cd4", paste0(cd4_subset, "_KEGG_up.xlsx"))
  )
}

if (!is.null(KEGG_down_df)) {
  write_xlsx(
    KEGG_down_df,
    file.path("GoKegg_cd4", paste0(cd4_subset, "_KEGG_down.xlsx"))
  )
}


## ---------------------------------------------------------------------------
## 6. Bubble plot helper (using utils plot_theme_common + save_plot)
## ---------------------------------------------------------------------------

plot_bubble <- function(df, title_text, filename) {
  if (is.null(df) || nrow(df) == 0) {
    return(invisible(NULL))
  }

  df_plot <- df[1:min(10, nrow(df)), ]
  df_plot$Description <- factor(df_plot$Description,
                                levels = rev(df_plot$Description))

  p <- ggplot(df_plot, aes(x = 1, y = Description)) +
    geom_point(aes(size = Count, color = -log10(pvalue))) +
    scale_color_gradient(low = "#56B1F7", high = "#CA0020") +
    scale_size(range = c(3, 8)) +
    labs(title = title_text, x = "", y = "") +
    plot_theme_common()

  save_plot(
    p,
    filename = file.path("GoKegg_cd4", filename),
    width    = 6,
    height   = 5
  )
}


## ---------------------------------------------------------------------------
## 7. Bubble plots for GO and KEGG (Up / Down)
## ---------------------------------------------------------------------------

plot_bubble(
  GO_up_df,
  title_text = paste0("GO Up: ", cd4_subset),
  filename   = paste0(cd4_subset, "_GO_Up_bubble.pdf")
)

plot_bubble(
  GO_down_df,
  title_text = paste0("GO Down: ", cd4_subset),
  filename   = paste0(cd4_subset, "_GO_Down_bubble.pdf")
)

plot_bubble(
  KEGG_up_df,
  title_text = paste0("KEGG Up: ", cd4_subset),
  filename   = paste0(cd4_subset, "_KEGG_Up_bubble.pdf")
)

plot_bubble(
  KEGG_down_df,
  title_text = paste0("KEGG Down: ", cd4_subset),
  filename   = paste0(cd4_subset, "_KEGG_Down_bubble.pdf")
)

message("CD4 GO/KEGG enrichment finished for subset: ", cd4_subset)

###############################################################################
# End of enrich_scRNA_cd4.R
###############################################################################
