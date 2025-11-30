# iCDC

# Project Structure

This repository contains the complete analysis pipeline for both **single-cell RNA-seq (scRNA-seq)** and **bulk RNA-seq** components of the study.  

---

## 📁 Repository Structure

```
main/
├── README.md                   # This file
│
├── bulk/
│   ├── README_bulk.md          # Bulk RNA-seq documentation
│   ├── bulk_1-2.R              # For T cell / DC bulk DEG enrichment
│   ├── bulk_3.R                # For whole-heart bulk DEG enrichment
│   ├── utils_bulk.R            # Shared bulk functions
│   ├── input/                  # (Empty)
│   └── output/                 # (Empty)
│
├── scRNA/
│   ├── README_scRNA.md         # scRNA-seq documentation
│   │
│   ├── utils/
│   │   └── utils_scRNA.R       # Shared functions for all scRNA workflows
│   │
│   ├── plots/                  # (Empty)
│   ├── results/                # (Empty)
│   │
│   ├── Main/
│   │   ├── whole_heart_main.R  # Whole-heart scRNA-seq pipeline
│   │   ├── cd45_main.R         # CD45+ immune compartments
│   │   └── cd4_T_main.R        # CD4+ T cells
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
│       ├── enrich_scRNA.R      # General IR vs Sham/DC unified enrichment
│       └── enrich_scRNA_cd4.R  # CD4-specific MI vs iCDC enrichment
│
└── LICENSE
```

---

## 🔒 Data Availability

To protect unpublished data, **all data and result folders are intentionally left empty**:

- `bulk/input/`
- `bulk/output/`
- `scRNA/plots/`
- `scRNA/results/`

Only the **full codebase** is provided to ensure complete reproducibility *without exposing sensitive information*.

---

## 📘 Notes

- All analysis scripts are written in **R**, validated end-to-end, and modularized via `utils_bulk.R` and `utils_scRNA.R`.
- The repository follows a strict structure to support future expansion and automated workflows.
- Subcluster modules, enrichment modules, and whole-heart/CD45/CD4 pipelines are fully separated for clarity.

---
