# Install the package versions used for the paper (R 4.4.1, Bioconductor
# 3.19; full record in sessionInfo.txt). Run once:
#   Rscript install_dependencies.R

options(repos = c(CRAN = "https://cloud.r-project.org"))
for (pkg in c("remotes", "BiocManager")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
BiocManager::install(version = "3.19", ask = FALSE, update = FALSE)
options(repos = BiocManager::repositories())

cran <- c(
  SeuratObject = "5.0.2",
  Seurat = "5.1.0",
  ggplot2 = "3.5.1",
  ggpubr = "0.6.0",
  dplyr = "1.1.4",
  circlize = "0.4.16",
  reshape2 = "1.4.4",
  RColorBrewer = "1.1-3",
  NMF = "0.27",
  ggalluvial = "0.12.5"
)
for (pkg in names(cran)) {
  if (!requireNamespace(pkg, quietly = TRUE) || packageVersion(pkg) != cran[[pkg]]) {
    remotes::install_version(pkg, version = cran[[pkg]], upgrade = "never")
  }
}

BiocManager::install(c("ComplexHeatmap", "BiocNeighbors"), ask = FALSE, update = FALSE)
remotes::install_github("immunogenomics/presto", upgrade = "never")
remotes::install_github("satijalab/azimuth@v0.5.0", upgrade = "never")
remotes::install_github("jinworks/CellChat@v2.1.2", upgrade = "never")
