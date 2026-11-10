#' Project paths
#'
#' Large inputs and cached intermediates are not tracked in git. By default
#' they are read from `01-InputData/` and `02-DataAnalysis/` under the project
#' root; both locations, and the external datasets, can be overridden with
#' environment variables (see `.Renviron.example`).

#' @describeIn input_path Path inside the input-data directory
#'   (`TOLDC_INPUT_DIR`, default `01-InputData`).
input_path <- function(...) file.path(Sys.getenv("TOLDC_INPUT_DIR", "01-InputData"), ...)

#' @describeIn input_path Path inside the cached-intermediates directory
#'   (`TOLDC_ANALYSIS_DIR`, default `02-DataAnalysis`).
analysis_path <- function(...) file.path(Sys.getenv("TOLDC_ANALYSIS_DIR", "02-DataAnalysis"), ...)

#' @describeIn input_path Output path for a figure panel (created on demand).
figure_path <- function(file) {
  dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
  file.path("results/figures", file)
}

#' @describeIn input_path Output path for a table (created on demand).
table_path <- function(file) {
  dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
  file.path("results/tables", file)
}

#' Resolve a path to an external resource that lives outside the project
#'
#' @param key One of `"biomart_gene_map"` (Ensembl -> symbol map used after
#'   Azimuth annotation), `"lin2020_dir"` (GSE154778 raw 10x output, one
#'   sub-directory per patient) or `"lin2020_meta"` (GSE154778 sample sheet).
#' @return The expanded path; errors if it is unset or missing.
external_path <- function(key) {
  env_var <- c(
    biomart_gene_map = "TOLDC_BIOMART_GENEMAP_PATH",
    lin2020_dir = "TOLDC_LIN2020_DIR",
    lin2020_meta = "TOLDC_LIN2020_META"
  )[[key]]
  path <- path.expand(Sys.getenv(env_var))
  if (!nzchar(path) || !file.exists(path)) {
    stop(sprintf("Set %s to the location of '%s' (see .Renviron.example). Current value: '%s'",
                 env_var, key, path))
  }
  path
}
