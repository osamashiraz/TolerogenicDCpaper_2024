#!/usr/bin/env Rscript
# Run every analysis script in order, each in its own R session.
# Arguments (e.g. --recompute) are passed through to every script.

scripts <- c(
  "scripts/01_azimuth_reference.R",
  "scripts/02_figure1.R",
  "scripts/03_figure3_cellchat.R",
  "scripts/04_external_validation.R"
)
args <- commandArgs(trailingOnly = TRUE)

for (script in scripts) {
  message("\n== ", script)
  status <- system2("Rscript", c(script, args))
  if (status != 0) stop(script, " failed (exit status ", status, ")")
}
