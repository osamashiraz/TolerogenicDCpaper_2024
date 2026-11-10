#' Differentially expressed genes of the CD141+ mDC subclusters
#'
#' Three comparisons (Methods a-c; Suppl. Tables S2-S4) share one
#' significance filter. Their intersection gives the Fig. 1D markers.

#' Significance filter for marker tables
#'
#' avg_log2FC >= 1, p_val_adj <= 0.05, pct.1 > 0.7 and pct.1 - pct.2 > 0.5.
#'
#' @param marker_df `FindMarkers()`/`FindAllMarkers()` output. A `diff`
#'   column (`pct.1 - pct.2`) is added if absent.
sig_filter <- function(marker_df, log2fc_min = 1, padj_max = 0.05,
                       pct1_min = 0.7, pct_diff_min = 0.5) {
  if (!"diff" %in% names(marker_df)) {
    marker_df$diff <- marker_df$pct.1 - marker_df$pct.2
  }
  subset(marker_df, avg_log2FC >= log2fc_min & p_val_adj <= padj_max &
           diff > pct_diff_min & pct.1 > pct1_min)
}

#' (a) Each mDC subcluster vs. the other two
#'
#' @param dc Output of [assign_dc_subset_labels()].
markers_subclusters_vs_each_other <- function(dc) {
  markers <- Seurat::FindAllMarkers(dc)
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$significant_DEG <- markers$gene %in% sig_filter(markers)$gene
  markers
}

#' (b) Each mDC subcluster vs. all other PDAC spleen cells
#'
#' @param pda PDAC object with subcluster labels in `cellTypes`
#'   ([merge_dc_subset_labels()]).
markers_subclusters_vs_all <- function(pda) {
  Seurat::Idents(pda) <- as.factor(pda$cellTypes)
  markers <- Seurat::FindAllMarkers(pda)
  markers <- subset(markers, cluster %in% mdc_subsets)
  markers$cluster <- as.character(markers$cluster)
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$significant_DEG <- markers$gene %in% sig_filter(markers)$gene
  markers
}

#' (c) The whole CD141+ mDC cluster vs. all other PDAC spleen cells
#'
#' Genes that separate the subclusters from each other (comparison a) are
#' not counted as significant here: this comparison is for shared markers.
#'
#' @param pda PDAC object with `predicted.cell_ontology_class`.
#' @param subcluster_sig Significant rows of comparison (a).
markers_whole_cluster_vs_rest <- function(pda, subcluster_sig) {
  markers <- Seurat::FindMarkers(
    pda, ident.1 = "cd141-positive myeloid dendritic cell", group.by = "predicted.cell_ontology_class"
  )
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$gene <- rownames(markers)
  sig <- subset(sig_filter(markers), gene %notin% subcluster_sig$gene)
  markers$significant_DEG <- markers$gene %in% sig$gene
  markers
}

#' Fig. 1D marker table (Suppl. Table S5)
#'
#' Subcluster-specific markers are genes significant in (a) and among the
#' top 50 per subcluster in (b); the 5 most significant whole-cluster
#' markers from (c) are labeled "shared". The top `top_n` genes per group
#' are kept (ties included, so a group can have more).
#'
#' @param each_other_sig,vs_all_sig,whole_sig `sig_filter()`-ed outputs of
#'   comparisons (a), (b) and (c).
build_figure1d_markers <- function(each_other_sig, vs_all_sig, whole_sig, n_shared = 5, top_n = 5) {
  # Many adjusted p-values underflow to 0, so row order decides ties: genes
  # are ranked by pct.1 - pct.2, and a gene significant in several
  # subclusters takes the row with the largest difference.
  each_other_sig <- each_other_sig[order(each_other_sig$diff, decreasing = TRUE), ]
  vs_all_sig <- vs_all_sig[order(vs_all_sig$diff, decreasing = TRUE), ]
  whole_sig <- whole_sig[order(whole_sig$p_val_adj), ]

  vs_all_top <- vs_all_sig %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = -p_val_adj, n = 50)
  shared <- head(whole_sig$gene, n_shared)
  genes <- unique(c(intersect(each_other_sig$gene, vs_all_top$gene), shared))

  res <- vs_all_sig[match(genes, vs_all_sig$gene), ]
  res <- res[order(res$p_val_adj, decreasing = TRUE), ]
  res[res$gene %in% shared, "cluster"] <- "shared"

  res <- res %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = -p_val_adj, n = top_n) %>%
    as.data.frame()
  res[, c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "diff", "gene")]
}
