# Shared setup for the analysis scripts. Run them from the project root:
#   Rscript scripts/02_figure1.R              # load cached intermediates
#   Rscript scripts/02_figure1.R --recompute  # rebuild them from the inputs

if (!file.exists("scripts/_setup.R")) {
  stop("Run the scripts from the project root, e.g. `Rscript scripts/run_all.R`.")
}
if (file.exists(".Renviron.local")) readRenviron(".Renviron.local")

options(toldc.recompute = "--recompute" %in% commandArgs(trailingOnly = TRUE) ||
          Sys.getenv("TOLDC_RECOMPUTE") == "1")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
})

# Palettes and labels first: other modules use them at load time.
source("scripts/R/colors_themes.R")
for (f in setdiff(list.files("scripts/R", pattern = "\\.R$", full.names = TRUE), "scripts/R/colors_themes.R")) {
  source(f)
}
