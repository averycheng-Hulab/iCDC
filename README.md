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

## Figure reproduction (Figs)

In addition to the demo pipeline, this repository provides **figure-level reproduction scripts** under `Figs/`.
These scripts use **minimal, de-identified figure inputs** (e.g., exported metadata tables, per-cell summary tables, gene expression long tables) to reproduce representative manuscript main figure panels without requiring full Seurat objects or large raw datasets.

### Figs: inputs and outputs

* **Figure inputs:** `Figs/data/`
  Small CSV tables exported from upstream analysis (Seurat/bulk results) for figure reproduction.
* **Figure scripts:** `Figs/scripts/`
  R scripts to reproduce figure panels (e.g., Fig3j, Fig3k).
* **Figure outputs:** `Figs/outputs/`
  Generated PDFs after running figure scripts.

### Run figure scripts

Run from the repository root:

```bash
Rscript Figs/scripts/Fig3j.R
Rscript Figs/scripts/Fig3k.R
```

Outputs will be saved to:

* `Figs/outputs/`

> Note: Figure scripts are designed to be lightweight and reproducible using the exported input tables.
> Full-resolution data and complete intermediate objects are not distributed in this GitHub repository.
> Reproduce (`Figs/`): figure scripts reproduce the manuscript results; differences from the manuscript are cosmetic (e.g., fonts, spacing, legend placement, PDF rendering). Final panel assembly was done in Adobe Illustrator.


---

## System requirements

This repository is a **collection of analysis scripts** (not a compiled software package).
To run the demos and/or reuse the analysis pipeline, you need:

* **R** (version to be specified by the user; see "Record your environment" below)
* Operating system: Linux / macOS / Windows
* Internet access is required for **first-time package installation** (CRAN/Bioconductor/GitHub)

### Dependencies (demo)

Installed by `install_packages_demo.R`:

* CRAN: `Seurat`, `ggplot2`, `dplyr`, `stringr`, `writexl`, `patchwork`
* Bioconductor: `clusterProfiler`, `org.Mm.eg.db`
* GitHub: `SCP` (`zhanghao-njmu/SCP`)

### Dependencies (full analysis)

Installed by `install_packages.R` (superset of demo):

* CRAN: `Seurat`, `harmony`, `ggplot2`, `patchwork`, `dplyr`, `stringr`, `clustree`, `writexl`
* Bioconductor: `clusterProfiler`, `org.Mm.eg.db`, `scRepertoire`, `enrichplot`
* GitHub: `DoubletFinder` (`chris-mcginnis-ucsf/DoubletFinder`), `SCP` (`zhanghao-njmu/SCP`)

---

## Repository structure

```
.
├── install_packages_demo.R        # install packages for demo
├── install_packages.R             # install packages for full analysis
├── run_demo.R                     # entry point to run bulk/scRNA demos
├── README.md                      # top-level documentation (this file)
├── LICENSE
├── CITATION.cff
├── .gitignore
├── demo/                          # lightweight demo inputs + generated outputs
│   ├── demo_bulk.R
│   ├── demo_scRNA.R
│   ├── bulk_input/
│   ├── bulk_output/               # (generated after running)
│   ├── scRNA_input/
│   └── scRNA_output/              # (generated after running)
├── Figs/                          # figure-level reproduction (lightweight)
│   ├── data/                      # exported figure inputs (CSV)
│   ├── scripts/                   # figure scripts (R)
│   └── outputs/                   # (generated) figure PDFs
├── bulk/                          # bulk RNA-seq analysis code
│   ├── DESeq2.R
│   ├── bulk_1-2.R
│   ├── bulk_3.R
│   ├── utils/
│   │   └── utils_bulk.R
│   └── README.md
└── scRNA/                         # single-cell RNA-seq analysis code
    ├── utils/
    │   └── utils_scRNA.R
    ├── ECM_genelist/
    │   └── ECMregulators.csv
    ├── Main/
    ├── Subcluster/
    ├── Enrichment/
    └── README.md
```

**Notes**

* `demo/*_input/` directories contain **lightweight demo inputs** so the pipeline can be run end-to-end.
* `demo/*_output/` directories are **generated after running** the demos.
* `Figs/data/` contains **minimal figure inputs** exported from upstream analysis for figure reproduction.
* Some directories may include `.gitkeep` files to preserve folder structure in GitHub.

---

## Overview of analyses

This repository implements two major transcriptomics analysis modules:

### Single-cell RNA-seq (scRNA)

Located under `scRNA/`.

Includes:

* Whole-heart multi-sample integration
* Immune-cell extraction (CD45)
* CD4 T-cell profiling
* Fibroblast re-clustering & ECM regulator program analysis
* Immune subclustering (Mac/Mono/DC, T/NK, Neutrophils, B-cells)
* CD4 T-cell TCR clonotype integration
* Unified GO/KEGG enrichment modules

All scripts share a unified analytical framework provided through `utils_scRNA.R`.

Detailed documentation is provided in:

```
scRNA/README.md
```

### Bulk RNA-seq

Located under `bulk/`.

In vitro T cell/DC and whole-heart bulk RNA-seq.

Detailed documentation is provided in:

```
bulk/README.md
```

---

## Data availability & privacy

This repository contains the **analysis code** together with **lightweight demo inputs** under `demo/` and **figure-level lightweight inputs** under `Figs/data/` to allow reviewers to run representative workflows and reproduce selected figure panels.

To protect unpublished and/or sensitive data and avoid substantial compute resources and time, full-resolution datasets and complete result folders are not distributed in this GitHub repository. In particular, folders such as:

* `demo/bulk_output/`
* `demo/scRNA_output/`
* `Figs/outputs/`

may be empty prior to running.

**Full data access information is provided in the manuscript “Data availability” section.**

---

## Software requirements

* **R ≥ 4.1**
* **Seurat ≥ 4.3**
* **Harmony ≥ 1.2.3**
* **DoubletFinder ≥ 2.0.4**
* **clusterProfiler ≥ 4.14.6**
* **scRepertoire ≥ 2.5.0**
* **tximport ≥ 1.32.0**
* **SCP ≥ 0.5.6**
* **ggplot2 / dplyr / stringr / patchwork / writexl / enrichplot / org.Mm.eg.db**

Each subfolder README contains additional script-specific requirements.

---

## Maintainer

**Guo Cheng**

Department of Cardiology, The Second Affiliated Hospital, School of Medicine, Zhejiang University
Research Center for Life Science and Human Health, Binjiang Institute of Zhejiang University

```

---
