# iCDC
# Bulk RNA-seq module

This folder contains all scripts used for bulk RNA-seq downstream analysis in this study.  
We distinguish two types of bulk datasets:

- **`bulk_1-2`**: T cell bulk RNA-seq and dendritic cell (DC) bulk RNA-seq  
  - DEG tables provided by the sequencing platform (FPKM-based reports).
  - We used the platform-supplied DESeq2 differential expression results as input.
- **`bulk_3`**: Whole-heart bulk RNA-seq  
  - Raw FASTQ files were re-processed by us from scratch.
  - Read alignment, quantification and DESeq2 analysis were run locally.
  - Final DEG tables are TPM-based reports from this in-house pipeline.

Although the upstream pipelines differ slightly (platform vs in-house),  
**all downstream functional analyses (GO/KEGG) in this repository start from gene-level DEGs reported by DESeq2**.  

---

## 1. Folder structure

```text
bulk/
├── bulk_1-2.R          # T cell & DC bulk DEG enrichment (FPKM-based tables)
├── bulk_3.R            # Whole-heart bulk DEG enrichment (TPM-based tables)
├── utils_bulk.R        # Shared enrichment and plotting functions
├── README_bulk.md      # This file
├── input/
│   ├── bulk_1-2/       # DEG Excel files for bulk_1-2
│   └── bulk_3/         # DEG CSV files for bulk_3
└── output/
    ├── bulk_1-2/       # DEG lists + GO/KEGG tables (bulk_1-2)
    └── bulk_3/         # DEG lists + GO/KEGG tables (bulk_3)
