# iCDC

This repository contains all analysis scripts for both **single-cell RNA-seq (scRNA-seq)** and **bulk RNA-seq** analysis pipelines used in this project.  
It is designed for transparent and reproducible analysis, while ensuring strict protection of unpublished data.

---
## 📁 Repository Structure

```
main/
├── README.md                         # Top-level project documentation
├── LICENSE                           
├── CITATION.cff
├── .gitignore                       
│
├── scRNA/                            # Single-cell RNA-seq workflows
│   ├── README.md                     # scRNA documentation
│   │
│   ├── utils/
│   │   └── utils_scRNA.R             # Shared functions for all scRNA scripts
│   │
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
    ├── utils                 
    │   └── utils_scRNA.R            # Shared bulk enrichment utilities
    │
    ├── bulk_1-2.R                   # in vitro cell (T/DC) bulk enrichment
    ├── bulk_3.R                     # Whole-heart bulk enrichment
    │
    ├── demo_input/                        
    └── demo_output/                      


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
in vitro T cell/DC bulk RNA-seq

### **bulk_3**  
Whole-heart bulk RNA-seq

Detailed documentation is provided in:

```

bulk/README.md

```

---

#  **Data Availability & Privacy**

Because this repository reflects **unpublished data**, directories intentionally left empty include:

- `scRNA/plots/`
- `scRNA/results/`
- `bulk/demo_input/`
- `bulk/demo_output/`

> **To protect unpublished data, all data and result folders are intentionally left empty.  
> Only code is provided to ensure reproducibility without exposing sensitive information.**

---

#  **Software Requirements**

- **R ≥ 4.1**
- **Seurat ≥ 4.3**
- **Harmony ≥ 1.2.3**
- **DoubletFinder ≥ 2.0.4**
- **clusterProfiler ≥ 4.14.6**
- **scRepertoire ≥ 2.5.0**
- **tximport ≥ 1.32.0**

- **ggplot2 / dplyr / stringr / patchwork /writexl /org.Mm.eg.db**

Each subfolder README contains additional script-specific requirements.

---

# 👤 Maintainer

**Guo Cheng**  

Department of Cardiology，The Second Affiliated Hospital, School of Medicine, Zhejiang University

Research Center for Life Science and Human Health, Binjiang Institute of Zhejiang University


