#!/usr/bin/env Rscript
# install_packages.R
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install installers
install.packages(c("BiocManager", "remotes"), quiet = TRUE)

cran_pkgs <- c(
  "Seurat","harmony","ggplot2","patchwork",
  "dplyr","stringr","clustree","writexl"
)

bioc_pkgs <- c("clusterProfiler","org.Mm.eg.db","scRepertoire")

# CRAN
to_install_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install_cran) > 0) install.packages(to_install_cran, dependencies = TRUE)

# Bioconductor
to_install_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

# GitHub
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", dependencies = TRUE, upgrade = "never")
}
if (!requireNamespace("SCP", quietly = TRUE)) {
  remotes::install_github("zhanghao-njmu/SCP", dependencies = TRUE, upgrade = "never")
}

cat("Done. Run demo:\n",
    "  Rscript run_demo.R --mode bulk\n",
    "  Rscript run_demo.R --mode scrna\n",
    "  Rscript run_demo.R --mode all\n", sep = "")
