#' Cell-type labels, color palettes and ggplot themes

#' Azimuth (Tabula Sapiens) cell-ontology class -> short label used in the paper
#'
#' Two classes present in the reference (cd1c+ mDC, pDC) are not listed: they
#' are not called in the spleen dataset after rare-class filtering.
ontology_to_label <- c(
  "memory b cell" = "B_mem cell",
  "naive b cell" = "B_naive cell",
  "plasma cell" = "Plasma cell",
  "cd4-positive, alpha-beta memory t cell" = "CD4+ T_mem cell",
  "cd8-positive, alpha-beta memory t cell" = "CD8+ T_mem cell",
  "cd8-positive, alpha-beta t cell" = "CD8+ T cell",
  "naive thymus-derived cd4-positive, alpha-beta t cell" = "CD4+ T_naive cell",
  "naive thymus-derived cd8-positive, alpha-beta t cell" = "CD8+ T_naive cell",
  "regulatory t cell" = "T_reg cell",
  "mature nk t cell" = "NK-T_mature cell",
  "type i nk t cell" = "NK-T_type-i cell",
  "nk cell" = "NK cell",
  "innate lymphoid cell" = "Innate Lymphoid cell",
  "classical monocyte" = "Monocyte",
  "intermediate monocyte" = "Intermediate Monocyte",
  "cd141-positive myeloid dendritic cell" = "CD141+ mDC",
  "macrophage" = "Macrophage",
  "neutrophil" = "Neutrophil",
  "platelet" = "Platelet",
  "endothelial cell" = "Endothelial cell",
  "erythrocyte" = "Erythrocyte",
  "hematopoietic stem cell" = "Stem cell"
)

#' Colors for the raw Tabula Sapiens ontology classes (Suppl. Fig. S2)
ontology_cols <- c(
  "memory b cell" = "#616161", "naive b cell" = "#546E7A", "plasma cell" = "#BDBDBD",
  "cd4-positive, alpha-beta memory t cell" = "#CDDC39",
  "cd8-positive, alpha-beta memory t cell" = "#8BC34A",
  "cd8-positive, alpha-beta t cell" = "#4CAF50",
  "naive thymus-derived cd4-positive, alpha-beta t cell" = "#FFF59D",
  "naive thymus-derived cd8-positive, alpha-beta t cell" = "#C5E1A5",
  "regulatory t cell" = "#64FFDA",
  "mature nk t cell" = "#9C27B0", "type i nk t cell" = "#BA68C8", "nk cell" = "#CE93D8",
  "innate lymphoid cell" = "#00ABA9",
  "classical monocyte" = "#0277BD", "intermediate monocyte" = "#4FC3F7",
  "cd141-positive myeloid dendritic cell" = "#E57373",
  "cd1c-positive myeloid dendritic cell" = "#E53935",
  "plasmacytoid dendritic cell" = "#EF9A9A",
  "macrophage" = "#FA6800", "neutrophil" = "#00ABA9", "platelet" = "#6D8764",
  "endothelial cell" = "#825A2C", "erythrocyte" = "#A20025",
  "hematopoietic stem cell" = "#E3C800"
)

#' Colors for the paper's cell-type labels, with the CD141+ mDC cluster split
#' into its three subclusters (mDC_1 = mregDC). The unsplit `"CD141+ mDC"`
#' label is kept for the few normal-spleen DCs that were not subclustered.
celltype_cols <- c(
  "B_mem cell" = "#616161", "B_naive cell" = "#546E7A", "Plasma cell" = "#BDBDBD",
  "CD4+ T_mem cell" = "#CDDC39", "CD8+ T_mem cell" = "#8BC34A", "CD8+ T cell" = "#4CAF50",
  "CD4+ T_naive cell" = "#FFEE58", "CD8+ T_naive cell" = "#C5E1A5", "T_reg cell" = "#64FFDA",
  "NK-T_mature cell" = "#9C27B0", "NK-T_type-i cell" = "#BA68C8", "NK cell" = "#CE93D8",
  "Innate Lymphoid cell" = "#00ABA9",
  "Monocyte" = "#0277BD", "Intermediate Monocyte" = "#4FC3F7",
  "CD141+ mDC_1" = "#D32F2F", "CD141+ mDC_2" = "#EF9A9A", "CD141+ mDC_3" = "#FFCDD2",
  "CD141+ mDC" = "#D32F2F",
  "Macrophage" = "#FA6800", "Neutrophil" = "#00ABA9", "Platelet" = "#6D8764",
  "Endothelial cell" = "#825A2C", "Erythrocyte" = "#A20025", "Stem cell" = "#E3C800"
)

#' Colors for the two cohorts
condition_cols <- c("Normal" = "#4FC3F7", "PDAC" = "pink")

#' The three CD141+ mDC subclusters
mdc_subsets <- paste0("CD141+ mDC_", 1:3)

#' Display order of cell types in CellChat / violin panels
celltype_order <- c(
  "CD141+ mDC", mdc_subsets, "Macrophage", "Monocyte",
  "T_reg cell", "CD4+ T_mem cell", "CD4+ T_naive cell", "CD8+ T cell", "CD8+ T_mem cell",
  "CD8+ T_naive cell", "Innate Lymphoid cell", "NK cell", "NK-T_mature cell", "NK-T_type-i cell",
  "B_mem cell", "B_naive cell", "Plasma cell", "Erythrocyte", "Endothelial cell", "Stem cell"
)

#' Blue-to-red gradient used for expression and score feature plots
expression_gradient <- function() {
  ggplot2::scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 4, name = "RdYlBu")))
}

#' Base ggplot theme for bar plots
myTheme <- ggplot2::theme_linedraw() + ggplot2::theme(
  text = ggplot2::element_text(size = 15),
  axis.text = ggplot2::element_text(size = 12, face = "bold", color = "darkgray"),
  axis.title = ggplot2::element_text(size = 13, face = "bold", color = "darkgray"),
  legend.text = ggplot2::element_text(size = 12, face = "bold", color = "#757575"),
  plot.subtitle = ggplot2::element_text(size = 8, face = "bold"),
  legend.box.background = ggplot2::element_blank(),
  legend.box.margin = ggplot2::margin(0, 0, 0, 0),
  legend.title = ggplot2::element_text(size = 13, face = "bold", color = "#757575"),
  rect = ggplot2::element_rect(linetype = "dotted", colour = "grey", linewidth = 0.1),
  panel.grid = ggplot2::element_line(linetype = "dotted"),
  legend.background = ggplot2::element_rect(fill = "white", linewidth = 4, colour = "white"),
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_line(color = "#BDBDBD", linetype = "dotted")
)

#' Negated `%in%`
`%notin%` <- Negate(`%in%`)
