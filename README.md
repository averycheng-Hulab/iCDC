# iCDC

This repository contains all analysis scripts for both **single-cell RNA-seq (scRNA-seq)** and **bulk RNA-seq** analysis pipelines used in this project.  
It is designed for transparent and reproducible analysis, while ensuring strict protection of unpublished data.

---
## 📁 Repository Structure

```
main/
├── README.md                         # Top-level project documentation
├── LICENSE                           # Open-source license
│
├── scRNA/                            # Single-cell RNA-seq workflows
│   ├── README.md                     # scRNA documentation
│   │
│   ├── utils/
│   │   └── utils_scRNA.R             # Shared functions for all scRNA scripts
│   │
│   ├── plots/                        # (Empty)
│   ├── results/                      # (Empty)
│   │
│   ├── Main/
│   │   ├── whole_heart_main.R        # Whole-heart scRNA integration workflow
│   │   ├── cd45_main.R               # CD45+ immune compartment analysis
│   │   └── cd4_T_main.R              # CD4+ T-cell analysis
│   │
│   ├── Subcluster/
│   │   ├── whole_heart_fibroblast.R        # Fibroblast re-clustering
│   │   ├── whole_heart_fibroblast_ecm.R    # Fibroblast ECM regulator program
│   │   ├── cd45_MacMonoDc_subcluster.R     # Monocyte / Mac / DC subclustering
│   │   ├── cd45_tnk_subcluster.R           # T/NK subclustering
│   │   ├── cd45_neutrophil_subcluster.R    # Neutrophil subclustering
│   │   ├── cd45_b_subcluster.R             # B-cell subclustering
│   │   └── cd4_T_TCR_module.R              # CD4 TCR integration / clonotype module
│   │
│   └── Enrichment/
│       ├── enrich_scRNA.R                  # General IR vs Sham/DC enrichment
│       └── enrich_scRNA_cd4.R              # CD4-specific MI vs iCDC enrichment
│
└── bulk/                            # Bulk RNA-seq workflows
    ├── README.md                    # Bulk RNA-seq documentation
    │
    ├── utils_bulk.R                 # Shared bulk enrichment utilities
    │
    ├── bulk_1-2.R                   # T-cell & DC bulk enrichment (platform FPKM DEGs)
    ├── bulk_3.R                     # Whole-heart bulk enrichment (in-house DESeq2)
    │
    ├── input/                       # (Empty)
    └── output/                      # (Empty)


```

---

#  **Overview of Analyses**

This repository implements two major transcriptomics analysis modules:

---

##  **Single-Cell RNA-seq (scRNA)**

Located under `scRNA/`.

Includes:

- Whole-heart multi-sample integration
- Immune-cell extraction (CD45)
- CD4 T-cell profiling
- Fibroblast re-clustering & ECM regulator program analysis
- Immune subclustering (Mac/Mono/DC, T/NK, Neutrophils, B-cells)
- CD4 T-cell TCR clonotype integration
- Unified GO/KEGG enrichment modules

All scripts share a unified analytical framework provided through `utils_scRNA.R`.

Detailed documentation is provided in:

```

scRNA/README.md

```

---

##  **Bulk RNA-seq**

Located under `bulk/`.

Includes two distinct data sources:

### **bulk_1–2**  
(T cell/DC bulk RNA-seq)

- DEGs tables
- GO/KEGG enrichment  
- Intersection sets where relevant  

### **bulk_3**  
(Whole-heart bulk RNA-seq)
  
 DEG tables  
- GO/KEGG enrichment  
- Trend analysis across contrasts  

Detailed documentation is provided in:

```

bulk/README.md

```

---

#  **Data Availability & Privacy**

Because this repository reflects **unpublished data**, directories intentionally left empty include:

- `scRNA/plots/`
- `scRNA/results/`
- `bulk/input/`
- `bulk/output/`

> **To protect unpublished data, all data and result folders are intentionally left empty.  
> Only code is provided to ensure reproducibility without exposing sensitive information.**

---

#  **Software Requirements**

- **R ≥ 4.1**
- **Seurat ≥ 4.3**
- **harmony**
- **DoubletFinder**
- **clusterProfiler**
- **org.Mm.eg.db**
- **writexl**
- **ggplot2 / dplyr / stringr / patchwork**

Each subfolder README contains additional script-specific requirements.

---

# 👤 Maintainer

**Guo Cheng**  

Department of Cardiology，The Second Affiliated Hospital, School of Medicine, Zhejiang University

Research Center for Life Science and Human Health, Binjiang Institute of Zhejiang University


