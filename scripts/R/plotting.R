#' Shared plot builders

#' Save a ggplot / patchwork object, or a base-graphics drawing, as PDF
#'
#' @param plot A ggplot-like object, or a function that draws on the
#'   current device (for ComplexHeatmap and CellChat chord diagrams).
#' @param file File name inside `results/figures/`.
#' @param png Also write a PNG copy (used for the figures shown in the README).
save_pdf <- function(plot, file, width = 7, height = 6, png = FALSE) {
  draw <- function() if (is.function(plot)) plot() else print(plot)
  path <- figure_path(file)
  grDevices::pdf(path, width = width, height = height)
  draw()
  grDevices::dev.off()
  message("  wrote ", path)
  if (png) {
    png_path <- sub("\\.pdf$", ".png", path)
    grDevices::png(png_path, width = width, height = height, units = "in", res = 150)
    draw()
    grDevices::dev.off()
    message("  wrote ", png_path)
  }
  invisible(path)
}

#' Write a data.frame as tab-separated text in `results/tables/`
save_table <- function(df, file) {
  path <- table_path(file)
  utils::write.table(df, path, row.names = FALSE, quote = FALSE, sep = "\t")
  message("  wrote ", path)
  invisible(path)
}

#' UMAP colored by a metadata column
umap_plot <- function(obj, group_by, cols, reduction = "umap", pt_size = 0.2, alpha = 0.8, ...) {
  Seurat::DimPlot(
    obj, group.by = group_by, reduction = reduction, cols = cols,
    pt.size = pt_size, alpha = alpha, order = TRUE, na.value = NA, ...
  )
}

#' Feature plot with the paper's blue-to-red expression gradient
feature_umap <- function(obj, features, reduction = "umap", pt_size = 0.8, limits = NULL, ...) {
  gradient <- if (is.null(limits)) {
    expression_gradient()
  } else {
    ggplot2::scale_colour_gradientn(
      colours = rev(RColorBrewer::brewer.pal(n = 4, name = "RdYlBu")), limits = limits
    )
  }
  Seurat::FeaturePlot(
    obj, features = features, reduction = reduction, order = TRUE,
    alpha = 0.8, pt.size = pt_size, ...
  ) & gradient
}

#' Cell-type composition per sample (Fig. 1B, Suppl. Figs. S3/S7/S8)
#'
#' @param obj Seurat object.
#' @param sample_col,celltype_col Metadata columns for sample and cell type.
#' @param y `"Freq"` for percent of the sample, `"n"` for cell counts.
#' @param celltypes Optional subset of cell types to show (e.g. the mDCs).
composition_bar_plot <- function(obj, sample_col, celltype_col, y = "Freq", celltypes = NULL) {
  counts <- obj@meta.data %>%
    dplyr::group_by(.data[[sample_col]], .data[[celltype_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(Freq = 100 * n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(Sample = dplyr::all_of(sample_col), `Cell Type` = dplyr::all_of(celltype_col))
  if (!is.null(celltypes)) {
    counts <- dplyr::filter(counts, `Cell Type` %in% celltypes)
  }

  ggplot2::ggplot(counts, ggplot2::aes(x = Sample, y = .data[[y]], fill = `Cell Type`)) +
    ggplot2::geom_bar(stat = "identity", show.legend = !is.null(celltypes)) +
    ggplot2::scale_fill_manual(values = celltype_cols) +
    myTheme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
      plot.margin = ggplot2::margin(0.1, 0.25, 0.1, 0.5, "cm")
    ) +
    ggplot2::xlab(" ") +
    ggplot2::ylab(if (y == "Freq") "Frequency %" else "Cell Count")
}

#' Stacked violin plots of a ligand-receptor gene pair across cell types
#'
#' @param obj Seurat object with cell-type identities.
#' @param genes Ligand and receptor gene symbols.
#' @param idents Cell types to show (kept in the identity-level order).
#' @param ylim Y-axis (log expression) limits.
lr_violin <- function(obj, genes, idents, ylim = c(1, 5.5)) {
  idents <- intersect(levels(Seurat::Idents(obj)), idents)
  Seurat::VlnPlot(
    obj, idents = idents, features = genes, cols = celltype_cols, ncol = 1, add.noise = FALSE,
    sort = FALSE, same.y.lims = TRUE, log = TRUE, layer = "data", pt.size = 0.01, adjust = 0.5
  ) &
    ggplot2::xlab(" ") &
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90, size = 10, face = "bold", vjust = 0.5, hjust = 1, color = celltype_cols[idents]
    )) &
    ggplot2::scale_y_continuous(name = "Log Expression", limits = ylim)
}
