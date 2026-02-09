
# ============================================================
# DESeq2 DEG analysis + merge exp(normalized expression) for two groups
# Inputs:
#   count_group1, count_group2: count matrices/data.frames (genes x samples)
#   exp_group1, exp_group2: FPKM/TPM matrices/data.frames (genes x samples)
#   They must have the same number of rows and represent the same set of genes,
#   Row names must be identical and in the same order across all four objects.
#   
#   frame_id.csv contains Ensembl gene IDs and corresponding gene symbols.
# ============================================================


suppressPackageStartupMessages({
  library(DESeq2)
})

############################################################
# Comparison group1: CAR vs IR
############################################################

group1_name <- "IR"
group2_name <- "CAR"

# DESeq2 contrast direction:

contrast_vec <- c("area", group2_name, group1_name)
out_file <- paste0(group2_name, "_vs_", group1_name, ".csv")

count_threshold  <- 10
min_total_counts <- 50

frame_id <- read.csv("frame_id.csv",stringsAsFactors = FALSE)
rownames(frame_id) <- frame_id$geneid_ensemble
# -----------------------------
# Prepare data for DESeq2
# -----------------------------
min_group1_samples <- floor(ncol(count_group1) / 2) + 1
min_group2_samples <- floor(ncol(count_group2) / 2) + 1

gene_id <- rownames(count_group1)
exp_display <- data.frame(
  gene_symbol = frame_id[gene_id,]$geneid_symbol,
  exp_group1,
  exp_group2,
  check.names = FALSE
)


count_mat_raw <- as.matrix(data.frame(count_group1, count_group2, check.names = FALSE))
suppressWarnings(storage.mode(count_mat_raw) <- "numeric")

keep_g1 <- rowSums(count_mat_raw[, 1:ncol(count_group1), drop = FALSE] > count_threshold) >= min_group1_samples
keep_g2 <- rowSums(count_mat_raw[, (ncol(count_group1) + 1):ncol(count_mat_raw), drop = FALSE] > count_threshold) >= min_group2_samples
keep_idx <- which(keep_g1 | keep_g2)

gene_id_keep  <- gene_id[keep_idx]
count_mat_raw <- count_mat_raw[keep_idx, , drop = FALSE]
exp_display   <- exp_display[keep_idx, , drop = FALSE]

gene_row_id <- make.unique(gene_id_keep)
rownames(count_mat_raw) <- gene_row_id
rownames(exp_display)   <- gene_row_id

# Gene-level counts aggregated from Salmon/tximport can be non-integer (estimated counts).
count_mat <- round(count_mat_raw)
storage.mode(count_mat) <- "integer"
sample_names_1 <- colnames(count_group1)
sample_names_2 <- colnames(count_group2)

if (is.null(sample_names_1) || any(sample_names_1 == "")) sample_names_1 <- paste0(group1_name, "_", seq_len(ncol(count_group1)))
if (is.null(sample_names_2) || any(sample_names_2 == "")) sample_names_2 <- paste0(group2_name, "_", seq_len(ncol(count_group2)))
colnames(count_mat) <- c(sample_names_1, sample_names_2)
sample_group <- c(rep(group1_name, ncol(count_group1)), rep(group2_name, ncol(count_group2)))
coldata <- data.frame(area = factor(sample_group, levels = c(group1_name, group2_name)))
rownames(coldata) <- colnames(count_mat)

# -----------------------------
# Run DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ area)
dds <- dds[rowSums(counts(dds)) > min_total_counts, ]

dds <- DESeq(dds)
res <- results(dds, contrast = contrast_vec)

# -----------------------------
# Format results and merge normalized expression
# -----------------------------
res_df <- as.data.frame(res)
res_df$gene_row_id <- rownames(res_df)

exp_aligned <- exp_display[res_df$gene_row_id, , drop = FALSE]

# Mean expression (FPKM/TPM) in each group
group1_exp_mean <- rowMeans(exp_aligned[, 2:(1 + ncol(exp_group1)), drop = FALSE], na.rm = TRUE)
group2_exp_mean <- rowMeans(exp_aligned[, (2 + ncol(exp_group1)):(1 + ncol(exp_group1) + ncol(exp_group2)), drop = FALSE], na.rm = TRUE)


out_df <- cbind(
  data.frame(
    gene_symbol = exp_aligned$gene_symbol,
    gene_row_id = res_df$gene_row_id,
    group1_exp_mean = group1_exp_mean,
    group2_exp_mean = group2_exp_mean,
    check.names = FALSE
  ),
  exp_aligned[, -1, drop = FALSE],
  res_df[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), drop = FALSE]
)

write.csv(out_df, file = out_file, row.names = FALSE)

message("Done. Output written to: ", out_file)
message("Contrast direction: log2FoldChange = ", group2_name, " / ", group1_name)



############################################################
# Comparison group2: CAR vs Vector
############################################################

group1_name <- "Vector"
group2_name <- "CAR"

# DESeq2 contrast direction:

contrast_vec <- c("area", group2_name, group1_name)
out_file <- paste0(group2_name, "_vs_", group1_name, ".csv")

count_threshold  <- 10
min_total_counts <- 50

frame_id <- read.csv("frame_id.csv",stringsAsFactors = FALSE)
rownames(frame_id) <- frame_id$geneid_ensemble
# -----------------------------
# Prepare data for DESeq2
# -----------------------------
min_group1_samples <- floor(ncol(count_group1) / 2) + 1
min_group2_samples <- floor(ncol(count_group2) / 2) + 1

gene_id <- rownames(count_group1)
exp_display <- data.frame(
  gene_symbol = frame_id[gene_id,]$geneid_symbol,
  exp_group1,
  exp_group2,
  check.names = FALSE
)


count_mat_raw <- as.matrix(data.frame(count_group1, count_group2, check.names = FALSE))
suppressWarnings(storage.mode(count_mat_raw) <- "numeric")

keep_g1 <- rowSums(count_mat_raw[, 1:ncol(count_group1), drop = FALSE] > count_threshold) >= min_group1_samples
keep_g2 <- rowSums(count_mat_raw[, (ncol(count_group1) + 1):ncol(count_mat_raw), drop = FALSE] > count_threshold) >= min_group2_samples
keep_idx <- which(keep_g1 | keep_g2)

gene_id_keep  <- gene_id[keep_idx]
count_mat_raw <- count_mat_raw[keep_idx, , drop = FALSE]
exp_display   <- exp_display[keep_idx, , drop = FALSE]

gene_row_id <- make.unique(gene_id_keep)
rownames(count_mat_raw) <- gene_row_id
rownames(exp_display)   <- gene_row_id

# Gene-level counts aggregated from Salmon/tximport can be non-integer (estimated counts).
count_mat <- round(count_mat_raw)
storage.mode(count_mat) <- "integer"
sample_names_1 <- colnames(count_group1)
sample_names_2 <- colnames(count_group2)

if (is.null(sample_names_1) || any(sample_names_1 == "")) sample_names_1 <- paste0(group1_name, "_", seq_len(ncol(count_group1)))
if (is.null(sample_names_2) || any(sample_names_2 == "")) sample_names_2 <- paste0(group2_name, "_", seq_len(ncol(count_group2)))
colnames(count_mat) <- c(sample_names_1, sample_names_2)
sample_group <- c(rep(group1_name, ncol(count_group1)), rep(group2_name, ncol(count_group2)))
coldata <- data.frame(area = factor(sample_group, levels = c(group1_name, group2_name)))
rownames(coldata) <- colnames(count_mat)

# -----------------------------
# Run DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ area)
dds <- dds[rowSums(counts(dds)) > min_total_counts, ]

dds <- DESeq(dds)
res <- results(dds, contrast = contrast_vec)

# -----------------------------
# Format results and merge normalized expression
# -----------------------------
res_df <- as.data.frame(res)
res_df$gene_row_id <- rownames(res_df)

exp_aligned <- exp_display[res_df$gene_row_id, , drop = FALSE]

# Mean expression (FPKM/TPM) in each group
group1_exp_mean <- rowMeans(exp_aligned[, 2:(1 + ncol(exp_group1)), drop = FALSE], na.rm = TRUE)
group2_exp_mean <- rowMeans(exp_aligned[, (2 + ncol(exp_group1)):(1 + ncol(exp_group1) + ncol(exp_group2)), drop = FALSE], na.rm = TRUE)


out_df <- cbind(
  data.frame(
    gene_symbol = exp_aligned$gene_symbol,
    gene_row_id = res_df$gene_row_id,
    group1_exp_mean = group1_exp_mean,
    group2_exp_mean = group2_exp_mean,
    check.names = FALSE
  ),
  exp_aligned[, -1, drop = FALSE],
  res_df[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), drop = FALSE]
)

write.csv(out_df, file = out_file, row.names = FALSE)

message("Done. Output written to: ", out_file)
message("Contrast direction: log2FoldChange = ", group2_name, " / ", group1_name)



