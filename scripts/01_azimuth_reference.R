# Tabula Sapiens spleen reference for Azimuth (Methods; Suppl. Fig. S2)
#
# Builds the Azimuth reference used to annotate the spleens (only with
# --recompute, or if it is missing) and plots the reference dataset.

source("scripts/_setup.R")

get_azimuth_reference()

ts_spleen <- readRDS(input_path("References", "TS_spleen", "cbab13ec-d454-4d3b-8ab0-86715b140d83.rds"))
save_pdf(
  umap_plot(ts_spleen, "cell_ontology_class", ontology_cols, pt_size = 0.8, alpha = 1),
  "SupplFigS2_TabulaSapiens_spleen_umap.pdf", width = 11, height = 6
)
