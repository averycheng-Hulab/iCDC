# iCDC
# Bulk RNA-seq module

This folder contains all scripts used for bulk RNA-seq downstream analysis in this study.  
We distinguish two types of bulk datasets:

- **`bulk_1-2`**: T cell bulk RNA-seq and dendritic cell (DC) bulk RNA-seq  
  - DEG tables (FPKM-based) provided by the sequencing platform.
- **`bulk_3`**: Whole-heart bulk RNA-seq 
  - DEG tables (TPM-based) reports from in-house pipeline.

Although the upstream pipelines differ slightly (platform vs in-house),  
**all downstream functional analyses (GO/KEGG) in this repository start from gene-level DEGs reported by DESeq2**.  

---

## 1. Folder structure

```text
bulk/
├── bulk_1-2.R          # T cell & DC bulk DEG enrichment
├── bulk_3.R            # Whole-heart bulk DEG enrichment
├── utils_bulk.R        # Shared enrichment and plotting functions
├── README_bulk.md      # This file
├── input/
│   ├── bulk_1-2/       # DEG tables bulk_1-2
│   └── bulk_3/         # DEG tables for bulk_3
└── output/
    ├── bulk_1-2/       # DEG lists + GO/KEGG tables + plots (bulk_1-2)
    └── bulk_3/         # DEG lists + GO/KEGG tables + plots (bulk_3)
