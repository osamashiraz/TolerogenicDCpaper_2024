#' PDAC subset and CD141+ mDC subclustering

#' Subset to PDAC spleens and re-process
#'
#' @param spleen_obj Annotated object from [label_cell_types()].
subset_pdac <- function(spleen_obj) {
  pda <- subset(spleen_obj, sample %in% "PDAC")
  pda <- Seurat::NormalizeData(pda)
  pda <- Seurat::FindVariableFeatures(pda)
  pda <- Seurat::ScaleData(pda, features = rownames(pda))
  pda <- Seurat::RunPCA(pda)
  pda <- Seurat::FindNeighbors(pda, dims = 1:15)
  pda <- Seurat::FindClusters(pda)
  Seurat::RunUMAP(pda, dims = 1:15)
}

#' Pick the number of PCs by the elbow heuristic
#'
#' First PC at which the cumulative share of standard deviation exceeds
#' `cum_var_threshold`% while that PC alone contributes under
#' `per_pc_threshold`% (hbctraining "elbow plot metric").
#'
#' @param stdev Per-PC standard deviations, e.g. `obj[["pca"]]@stdev`.
select_elbow_pcs <- function(stdev, cum_var_threshold = 90, per_pc_threshold = 5) {
  pct <- stdev / sum(stdev) * 100
  which(cumsum(pct) > cum_var_threshold & pct < per_pc_threshold)[1]
}

#' Re-cluster the CD141+ mDC population of the PDAC spleens
#'
#' @param pda Output of [subset_pdac()].
#' @return The mDC-only object, clustered at resolutions 0.1-0.8.
subcluster_mdc <- function(pda) {
  dc <- subset(pda, cellTypes %in% "CD141+ mDC")
  dc <- Seurat::NormalizeData(dc)
  dc <- Seurat::FindVariableFeatures(dc, selection.method = "dispersion")
  dc <- Seurat::ScaleData(dc, features = rownames(dc))
  dc <- Seurat::RunPCA(dc)

  n_pcs <- select_elbow_pcs(dc[["pca"]]@stdev)
  dc <- Seurat::FindNeighbors(dc, dims = 1:n_pcs)
  dc <- Seurat::FindClusters(dc, resolution = seq(0.1, 0.8, 0.1))
  dc <- Seurat::RunUMAP(dc, dims = 1:n_pcs)
  Seurat::UpdateSeuratObject(dc)
}

#' Name the mDC subclusters
#'
#' At resolution 0.2 the published clustering has three clusters; raw
#' cluster 2 is the IDO1-high mregDC population (mDC_1). Cluster IDs are
#' arbitrary, so this mapping must be re-checked if the subclustering is
#' recomputed (e.g. with `IDO1` expression).
#'
#' @param dc Output of [subcluster_mdc()].
#' @param mapping Raw cluster ID -> subcluster number; unlisted IDs map to
#'   `.default`.
#' @return `dc` with a `dc_subsets` column, also set as the identity.
assign_dc_subset_labels <- function(dc, mapping = c("2" = "1", "1" = "2", .default = "3")) {
  raw <- as.character(dc$RNA_snn_res.0.2)
  suffix <- ifelse(raw %in% names(mapping), mapping[raw], mapping[[".default"]])
  dc$dc_subsets <- paste0("CD141+ mDC_", suffix)
  Seurat::Idents(dc) <- factor(dc$dc_subsets, levels = mdc_subsets)
  dc
}

#' Replace the single CD141+ mDC label with the subcluster labels
#'
#' @param obj Any object containing the subclustered cells (PDAC or all
#'   spleens).
#' @param dc Output of [assign_dc_subset_labels()].
merge_dc_subset_labels <- function(obj, dc) {
  cells <- intersect(colnames(obj), colnames(dc))
  obj$cellTypes[cells] <- dc$dc_subsets[cells]
  obj
}
