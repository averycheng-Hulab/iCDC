# iCDC

📌 Integrated scRNA-seq & Bulk RNA-seq Analysis Pipeline

This repository contains the full analysis framework used in this study, including:

Single-cell RNA-seq (scRNA-seq) processing

Immune subcluster workflows

Functional enrichment analysis (GO/KEGG)

Bulk RNA-seq DEG analysis

Unified utilities for reproducibility

The goal of this repository is to provide fully reproducible pipelines without exposing any unpublished biological data.

🔒 Data Availability & Confidentiality

To protect unpublished data, all data and result folders are intentionally left empty.
Only analysis code is provided to ensure full reproducibility without leaking sensitive biological information.

The following folders are intentionally empty:

scRNA/plots/

scRNA/results/

bulk/input/

bulk/output/

data_raw/ (if used)

data_processed/ (if used)

Users should populate these folders with their own data following the same structure.

📁 Repository Structure
main/
├── README.md                      # This file
│
├── bulk/
│   ├── README_bulk.md             # Bulk RNA-seq documentation
│   ├── bulk_1-2.R                 # For T cell / DC bulk DEG enrichment
│   ├── bulk_3.R                   # For whole-heart bulk DEG enrichment
│   ├── utils_bulk.R               # Shared bulk functions
│   ├── input/                     # (Empty; user must provide DEG files)
│   └── output/                    # (Empty; enrichment results)
│
├── scRNA/
│   ├── README_scRNA.md            # scRNA-seq documentation
│   │
│   ├── utils/
│   │   └── utils_scRNA.R          # Shared functions for all scRNA workflows
│   │
│   ├── plots/                     # (Empty; for storing figures)
│   ├── results/                   # (Empty; for storing analysis outputs)
│   │
│   ├── Main/
│   │   ├── whole_heart_main.R     # Whole-heart scRNA-seq pipeline
│   │   ├── cd45_main.R            # CD45+ immune compartments
│   │   └── cd4_T_main.R           # CD4+ T cells
│   │
│   ├── Subcluster/
│   │   ├── whole_heart_fibroblast.R
│   │   ├── whole_heart_fibroblast_ecm.R
│   │   ├── cd45_MacMonoDc_subcluster.R
│   │   ├── cd45_tnk_subcluster.R
│   │   ├── cd45_neutrophil_subcluster.R
│   │   └── cd45_b_subcluster.R
│   │
│   └── Enrichment/
│       ├── enrich_scRNA.R         # General IR vs Sham/DC unified enrichment
│       └── enrich_scRNA_cd4.R     # CD4-specific MI vs iCDC enrichment
│
└── LICENSE (optional)

🧬 Overview of Analysis Modules
1️⃣ scRNA-seq Main Pipelines

Located in scRNA/Main/, these scripts perform:

Quality control (QC)

Doublet removal

SCTransform normalization

PCA / Harmony batch correction

Clustering, UMAP/tSNE visualization

Cell type annotation (supervised rules + marker genes)

Object saving for downstream analysis

2️⃣ Subcluster Workflows

Located in scRNA/Subcluster/, handling:

Fibroblast & ECM modules

CD45⁺ immune compartments (T/NK, B, Neutrophil, MacMonoDC)

Reclustering and refined annotation

Marker-based subtype inference

All subcluster scripts rely on utils_scRNA.R for standardized processing.

3️⃣ Enrichment Analysis

Located in scRNA/Enrichment/:

enrich_scRNA.R: Generic IR vs Sham/DC enrichment across any celltype

enrich_scRNA_cd4.R: CD4-specific MI vs iCDC analysis

GO/KEGG enrichment performed via clusterProfiler

Bubble plots output to folder

4️⃣ Bulk RNA-seq Module

Located in bulk/:

bulk_1-2: T cell/DC bulk DEG (platform-generated)

bulk_3: Whole-heart bulk DEG (in-house pipeline)

All downstream GO/KEGG analyses standardized to DESeq2 DEG input

Functions centralized in utils_bulk.R

▶️ How to Run
scRNA-seq
source("scRNA/Main/whole_heart_main.R")
source("scRNA/Subcluster/whole_heart_fibroblast.R")
source("scRNA/Enrichment/enrich_scRNA.R")

Bulk RNA-seq
source("bulk/bulk_1-2.R")
source("bulk/bulk_3.R")

Enrichment example
run_enrich_scRNA(
    obj = fib_sl_object,
    celltype = "F-SL",
    outdir = "enrich_FSL"
)

🤝 Reproducibility Notes

All pipelines require R ≥ 4.2

All scripts use SCTransform and Harmony for batch correction

All cluster annotations follow fixed mapping rules (documented in subcluster scripts)

Code is modular and function-driven using utils_scRNA.R & utils_bulk.R
