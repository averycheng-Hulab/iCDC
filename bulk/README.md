# iCDC

---

# Bulk RNA-seq Analysis

This directory contains all scripts used for the bulk RNA-seq workflows in this project,  in vitro cell (T/DC) bulk datasets and whole-heart bulk RNA-seq analysis.

All analyses follow a unified GO/KEGG enrichment framework and use shared utilities for Mus musculus functional annotation.

### Bulk workflow (DEG sources)

- **bulk_1-2 (in vitro cells):** Differential expression results were obtained from the sequencing provider’s DESeq2 standardized reporting pipeline (downloaded DEG tables). This repository reproduces the downstream analyses starting from these DEG tables.
- **bulk_3 (heart tissue):** DEGs were generated in-house by running DESeq2 from count matrices, followed by the same downstream framework.

---

# 📁 Directory Structure

```

bulk/
├── README.md                     # This file
│
├── utils/
│   └── utils_bulk.R              # Shared enrichment + plotting utilities
│                     
├── bulk_1-2.R                    # in vitro cell (T/DC) bulk enrichment
└── bulk_3.R                      # Whole-heart bulk enrichment


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
source(file.path("bulk","utils","utils_bulk.R"))
````

---

# **2. Bulk Analysis Modules**

---

## **2.1 T cell / DC Bulk RNA-seq (bulk_1-2.R)**

This module processes  **DESeq2 outputs** . 

It performs:

* Load DEG tables
* DEG filtering
* GO enrichment
* KEGG enrichment
* Save:

  * Up/Down DEG lists
  * GO enrichment tables
  * KEGG enrichment tables
  * Example plots for each comparison

The FPKM filtering (expression cutoff) is **only used in bulk_1-2**.

---

## **2.2 Whole-Heart Bulk RNA-seq (bulk_3.R)**

This module processes **DESeq2 outputs**.

It performs:

* Load DEG tables
* DEG filtering
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



---

# **3. Input Requirements**

---

## **3.1 bulk_1-2: FPKM-based Excel tables**

Required columns:

| Column                 | Description      |
| ---------------------- | ---------------- |
| `gene_symbol`          | Gene symbol      |
| `log2FoldChange`       | Fold-change      |
| `padj`                 | Adjusted p-value |
| Two expression columns | For FPKM cutoff  |

You must update group-average expression column indices:

```r
expr_col_idx <- c(8, 9)
```

---

## **3.2 bulk_3: TPM-based CSV DESeq2 tables**

Required columns:

| Column           | Description      |
| ---------------- | ---------------- |
| `gene_symbol`    | Gene symbol      |
| `log2FoldChange` | DESeq2 FC        |
| `pvalue`         | Raw p-value      |
| `padj`           | Adjusted p-value |
| TPM columns      | One or more      |


---

# **4. Output Structure**

Each bulk module produces:

* DEG up/down lists
* GO enrichment tables
* KEGG enrichment tables
* Intersection DEG sets (bulk_3)
* Trend gene sets (bulk_3)
* Example plots for upregulated GO/KEGG

---

