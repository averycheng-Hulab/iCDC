# iCDC

---

# Bulk RNA-seq Analysis (Bulk)

This directory contains all scripts used for the bulk RNA-seq workflows in this project, including platform-derived T cell / DC bulk datasets and in-house whole-heart bulk RNA-seq analysis.

All analyses follow a unified GO/KEGG enrichment framework and use shared utilities for Mus musculus functional annotation.

---

# 📁 Directory Structure

```

bulk/
├── README_bulk.md                     # This file
│
├── utils_bulk.R                       # Shared enrichment + plotting utilities
│
├── bulk_1-2.R                         # T cell & DC bulk enrichment (FPKM-based)
├── bulk_3.R                           # Whole-heart bulk enrichment (TPM-based)
│
├── input/                             # (Empty)
│   ├── bulk_1-2/                      # (Empty) platform Excel DEG tables
│   └── bulk_3/                        # (Empty) in-house DESeq2 CSV DEG tables
│
└── output/                            # (Empty)
├── bulk_1-2/                      # (Empty) enrichment outputs
└── bulk_3/                        # (Empty) enrichment outputs

````

---

# **1. Utilities**

### **utils_bulk.R**

This file provides shared helper functions used by all bulk analysis scripts:

* GO Biological Process enrichment (clusterProfiler)
* KEGG pathway enrichment (mmu)
* Automatic ID conversion (SYMBOL ↔ ENTREZ)
* Reconstruction of KEGG gene lists in SYMBOL format
* Unified and error-safe enrichment wrappers
* Basic visualization modules:
  - Barplot for GO Biological Process
  - Dot/Bubble plot for KEGG pathways

All bulk workflows import this file:

```r
source(file.path("bulk", "utils_bulk.R"))
````

---

# **2. Bulk Analysis Modules**

---

## **2.1 T cell / DC Bulk RNA-seq (bulk_1-2.R)**

This module processes **platform-generated FPKM-based** DEG tables (Excel format). It performs:

* Load Excel DEG tables (platform DESeq2 results)
* Apply **FPKM expression cutoff** (platform-only requirement)
* DEG filtering:

  * Upregulated: log2FC > 1
  * Downregulated: log2FC < –1
* GO enrichment
* KEGG enrichment
* Save:

  * Up/Down DEG lists
  * GO enrichment tables
  * KEGG enrichment tables
  * Example plots for each comparison

The FPKM filtering (expression cutoff) is **only used in bulk_1-2**, and this is explained in this README.

Output directory:

```
bulk/output/bulk_1-2/
```

---

## **2.2 Whole-Heart Bulk RNA-seq (bulk_3.R)**

This module processes **in-house TPM-based DESeq2 outputs** generated from raw FASTQ files.

It performs:

* Load CSV DEG tables
* DEG selection criteria:

  * Upregulated: log2FC > 1 & pvalue < 0.01
  * Downregulated: log2FC < –1 & pvalue < 0.01
* GO enrichment
* KEGG enrichment
* Compute intersection and trend modules:

  * IR_vs_CAR ∩ Vector_vs_CAR (up/down)
  * Sham_vs_IR ∩ Sham_vs_Vector (up/down)
  * Trend:
    *Up in CAR-related & Down in Sham-related*
    *Down in CAR-related & Up in Sham-related*
* Produce:

  * DEG tables (up/down)
  * Enrichment tables
  * Trend gene sets + enrichment tables
  * RData files (for further downstream figures)

Output directory:

```
bulk/output/bulk_3/
```

---

# **3. Input Requirements**

---

## **3.1 bulk_1-2: Platform FPKM-based Excel tables**

Required columns:

| Column                 | Description      |
| ---------------------- | ---------------- |
| `gene_name`            | Gene symbol      |
| `log2FoldChange`       | Fold-change      |
| `padj`                 | Adjusted p-value |
| Two expression columns | For FPKM cutoff  |

You must update expression column indices:

```r
expr_col_idx <- c(8, 9)
```

---

## **3.2 bulk_3: In-house TPM-based CSV DESeq2 tables**

Required columns:

| Column           | Description      |
| ---------------- | ---------------- |
| `geneid_symbol`  | Gene symbol      |
| `log2FoldChange` | DESeq2 FC        |
| `pvalue`         | Raw p-value      |
| `padj`           | Adjusted p-value |
| TPM columns      | One or more      |

**No TPM expression cutoff is used**, as DESeq2 uses raw counts and statistical filtering.

---

# **4. Output Structure**

Each bulk module produces:

* DEG up/down lists
* GO enrichment tables
* KEGG enrichment tables
* Intersection DEG sets (bulk_3)
* Trend gene sets (bulk_3)
* Example plots for upregulated GO/KEGG

These are saved into the corresponding output folders:

```
bulk/output/bulk_1-2/
bulk/output/bulk_3/
```

These folders are intentionally empty in the repository.

---

# **5. Data Availability and Privacy**

Because this repository is associated with **unpublished work**, no data files are included:

* `bulk/input/`
* `bulk/output/`
* All platform Excel DEG tables
* All in-house DESeq2 CSV tables
* Any TPM/raw FASTQ dependencies

Only analysis **code** is provided.

> **To protect unpublished data, all data and result folders are intentionally left empty.
> Only code is provided to ensure reproducibility without exposing sensitive information.**

---

# **6. Software Requirements**

* R ≥ 4.1
* clusterProfiler
* org.Mm.eg.db
* enrichplot
* dplyr
* ggplot2
* readxl / writexl

---

# **7. Maintainer: Guo Cheng**

```

---

