# iCDC

This repository contains all analysis scripts for both **single-cell RNA-seq (scRNA-seq)** and **bulk RNA-seq** analysis pipelines used in this project.  
It is designed for transparent and reproducible analysis, while ensuring strict protection of unpublished data.


---

## Quick start (demo)

Please run demos **from the repository root**.

### 1) Install dependencies (demo only)

For running the **lightweight demos** in `demo/`, it is sufficient to install the demo dependencies:

```bash
Rscript install_packages_demo.R
````

> `install_packages.R` installs additional packages used in the **full analysis** (not required for the demo).
> If you only want to run the demo, you can skip `install_packages.R`.

### 2) Run the demo


Bulk demo:

```bash
Rscript run_demo.R --mode bulk
```

scRNA demo:

```bash
Rscript run_demo.R --mode scrna
```

Run both:

```bash
Rscript run_demo.R --mode all
```

### 3) Inputs and outputs

* **Bulk demo inputs:** `demo/bulk_input/`
  **Bulk demo outputs:** `demo/bulk_output/`

* **scRNA demo inputs:** `demo/scRNA_input/` (e.g., `seurat_demo.rds`)
  **scRNA demo outputs:** `demo/scRNA_output/`

> Note: The demo is intended to validate the pipeline and generate representative outputs.
> It is not designed to reproduce all manuscript figures due to compute/time constraints.



---
## Repository structure
---

```
.
├── install_packages_demo.R        # install packaged for demo
├── install_packages.R             # install packaged for full analysis
├── run_demo.R                     # entry point to run bulk/scRNA demos
├── README.md                      # top-level documentation (this file)
├── LICENSE
├── CITATION.cff
├── .gitignore
├── demo/                          # lightweight demo inputs + generated outputs
│   ├── demo_bulk.R                # bulk demo script (reads demo/bulk_input → writes demo/bulk_output)
│   ├── demo_scRNA.R               # scRNA demo script (reads demo/scRNA_input → writes demo/scRNA_output)
│   ├── bulk_input/                # demo DEG CSV inputs (small example files)
│   ├── bulk_output/               # demo bulk outputs (generated after running)
│   ├── scRNA_input/               # demo Seurat object(s), e.g., demo_fib.rds (downsampled)
│   └── scRNA_output/              # demo scRNA outputs (generated after running)
├── bulk/                          # bulk RNA-seq analysis code
│   ├── bulk_1-2.R
│   ├── bulk_3.R
│   ├── utils/
│   │   └── utils_bulk.R
│   ├── demo_input/                # (optional) alternative bulk demo input location
│   ├── demo_output/               # (optional) alternative bulk demo output location
│   └── README.md
└── scRNA/                         # single-cell RNA-seq analysis code
    ├── utils/
    │   └── utils_scRNA.R
    ├── ECM_genelist/
    │   └── ECMregulators.csv
    ├── Main/                      # main analysis scripts
    ├── Subcluster/                # subclustering scripts
    ├── Enrichment/                # enrichment scripts
    └── README.md
```

**Notes**

* The `demo/*_input/` directories contain **lightweight demo inputs** (small DEG tables and a downsampled Seurat object) so the pipeline can be run end-to-end.
* The `demo/*_output/` directories are **generated after running** the demos (placeholders may exist before execution).
* Some directories may include `.gitkeep` files to preserve folder structure in GitHub.

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

### **bulk_1–2**  
in vitro T cell/DC bulk RNA-seq

### **bulk_3**  
Whole-heart bulk RNA-seq

Detailed documentation is provided in:

```

bulk/README.md

```

---

# **Data Availability & Privacy**

This repository contains the **analysis code** together with **lightweight demo inputs** under `demo/` to allow reviewers to run the general workflow. The demo materials are **downsampled and de-identified** and are intended for demonstrating the pipeline rather than reproducing the full manuscript figures.

To protect unpublished and/or sensitive data and avoid substantial compute resources and time , full-resolution datasets and complete result folders are not distributed in this GitHub repository. In particular, folders such as:

- `demo/bulk_output/` (demo outputs are generated under here)
- `demo/scRNA_output/` ( demo outputs are generated under here)

may be empty prior to running the demos.

**Full data access information is provided in the manuscript “Data availability” section.**

---

#  **Software Requirements**

- **R ≥ 4.1**
- **Seurat ≥ 4.3**
- **Harmony ≥ 1.2.3**
- **DoubletFinder ≥ 2.0.4**
- **clusterProfiler ≥ 4.14.6**
- **scRepertoire ≥ 2.5.0**
- **tximport ≥ 1.32.0**
- **tximport ≥ 1.32.0**
- **SCP ≥ 0.5.6**

- **ggplot2 / dplyr / stringr / patchwork /writexl /org.Mm.eg.db**

Each subfolder README contains additional script-specific requirements.

---

# 👤 Maintainer

**Guo Cheng**  

Department of Cardiology，The Second Affiliated Hospital, School of Medicine, Zhejiang University

Research Center for Life Science and Human Health, Binjiang Institute of Zhejiang University


