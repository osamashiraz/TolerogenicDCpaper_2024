#' Tolerogenic ("mregDC") gene signature (Fig. 1E)

#' The 17-gene mregDC-associated signature defined in the Methods
toldc_geneset <- c(
  "IDO1", "IDO2", "CCR7", "CCL19", "CD274", "CD86", "CLEC2D", "CDH1", "AHR",
  "TXN", "HMOX1", "NECTIN2", "TGFB1", "CD163", "CXCL13", "IFNGR2", "NDUFA2"
)

#' Score a gene set per cell with `Seurat::AddModuleScore()`
#'
#' @param obj Seurat object.
#' @param genes Character vector of gene symbols.
#' @param name Name of the resulting metadata column.
#' @param ctrl Number of control genes per expression bin.
#' @return `obj` with a `name` column holding the module score.
add_signature_score <- function(obj, genes = toldc_geneset, name = "tol_DC", ctrl = 50) {
  obj <- Seurat::AddModuleScore(obj, features = list(genes), ctrl = ctrl, name = name, seed = 123)
  obj[[name]] <- obj[[paste0(name, "1")]]
  obj[[paste0(name, "1")]] <- NULL
  obj
}

#' Violin plot of the signature score with two-sided t-tests (Fig. 1E)
#'
#' Compares mDC_1 vs. mDC_2 and mDC_3 vs. mDC_2 (`ggpubr::stat_compare_means`,
#' `method = "t.test"`); comparisons with a group absent from `idents` are
#' skipped.
#'
#' @param obj Object with a `tol_DC` score and cell-type identities.
#' @param idents Identities to show.
signature_violin <- function(obj, idents = mdc_subsets) {
  idents <- intersect(levels(Seurat::Idents(obj)), idents)
  comparisons <- Filter(
    function(pair) all(pair %in% idents),
    list(c("CD141+ mDC_1", "CD141+ mDC_2"), c("CD141+ mDC_3", "CD141+ mDC_2"))
  )
  p <- Seurat::VlnPlot(
    obj, idents = idents, features = "tol_DC", cols = celltype_cols, add.noise = FALSE,
    sort = FALSE, pt.size = 0.01, adjust = 0.5
  ) & ggplot2::xlab("") & ggplot2::scale_y_continuous(name = "Signature\nScore")
  if (length(comparisons) > 0) {
    p <- p & ggpubr::stat_compare_means(comparisons = comparisons, label = "p.signif", method = "t.test")
  }
  p
}
