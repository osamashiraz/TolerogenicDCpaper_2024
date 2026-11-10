#' mDC subclusters and key ligand-receptor pairs in other datasets
#'
#' Each query dataset is mapped onto the in-house PDAC spleen object with
#' Seurat label transfer (`FindTransferAnchors()` / `TransferData()`), then
#' the key LR pairs and the tolerogenic signature are plotted per transferred
#' cell type (Fig. 1E right; Suppl. Figs. S7, S8).

#' Violin panels for the key LR pairs. Unlike [lr_specs] there is no CellChat
#' network here, so every cell list is explicit.
external_lr_specs <- local({
  checkpoint_targets <- c(mdc_subsets, "T_reg cell", "CD4+ T_mem cell", "CD4+ T_naive cell",
                          "CD8+ T cell", "CD8+ T_mem cell", "CD8+ T_naive cell", "NK cell",
                          "NK-T_mature cell", "NK-T_type-i cell", "Plasma cell", "Stem cell")
  specs <- Filter(function(s) s$name != "IDO1-AHR", lr_specs)
  lapply(specs, function(s) {
    if (identical(s$cells, "interacting")) s$cells <- checkpoint_targets
    s
  })
})

#' Load the Lin et al. 2020 PDAC dataset (GSE154778) and cluster it
#'
#' One 10x directory per patient, merged and clustered as a whole. The
#' immune compartment is later taken from clusters at resolution 0.1.
#'
#' @param data_dir,meta_path See `external_path("lin2020_dir")`,
#'   `external_path("lin2020_meta")`. The sample sheet has columns
#'   `Patient`, `ID`, `Type`, `Site`, `Age`, `Gender`, `Stage`, `Grade`.
load_lin2020 <- function(data_dir, meta_path) {
  meta <- utils::read.delim(meta_path)
  rownames(meta) <- meta$Patient

  samples <- lapply(meta$Patient, function(id) {
    obj <- Seurat::CreateSeuratObject(
      counts = Seurat::Read10X(file.path(data_dir, id)), min.cells = 3, min.features = 200, project = id
    )
    for (col in c("Patient", "ID", "Type", "Site", "Age", "Gender", "Stage", "Grade")) {
      obj[[if (col == "Patient") "patient" else col]] <- meta[id, col]
    }
    obj
  })
  obj <- merge(samples[[1]], y = samples[-1], add.cell.ids = meta$Patient,
               project = "2020_han_pda_met", merge.data = TRUE)
  obj <- SeuratObject::JoinLayers(obj)

  obj <- Seurat::NormalizeData(obj)
  obj <- Seurat::FindVariableFeatures(obj)
  obj <- Seurat::ScaleData(obj, features = rownames(obj))
  obj <- Seurat::RunPCA(obj)
  obj <- Seurat::FindNeighbors(obj, dims = 1:10)
  obj <- Seurat::FindClusters(obj, resolution = seq(0, 1, 0.1))
  Seurat::RunUMAP(obj, dims = 1:30)
}

#' Re-process a query dataset and transfer the in-house cell-type labels
#'
#' @param obj Query Seurat object.
#' @param reference PDAC spleen object with subcluster labels in `cellTypes`.
#' @param resolution Clustering resolution(s) (`NULL` = Seurat default).
#' @param scale_all Scale all genes (`TRUE`) or only variable genes.
#' @param min_cells Transferred cell types with fewer cells are dropped.
transfer_labels <- function(obj, reference, resolution = NULL, scale_all = TRUE, min_cells = 5) {
  obj <- Seurat::NormalizeData(obj)
  obj <- Seurat::FindVariableFeatures(obj)
  obj <- Seurat::ScaleData(obj, features = if (scale_all) rownames(obj) else NULL)
  obj <- Seurat::RunPCA(obj)
  obj <- Seurat::FindNeighbors(obj, dims = 1:15)
  obj <- if (is.null(resolution)) Seurat::FindClusters(obj) else Seurat::FindClusters(obj, resolution = resolution)
  obj <- Seurat::RunUMAP(obj, dims = 1:30)

  anchors <- Seurat::FindTransferAnchors(reference, obj)
  predictions <- Seurat::TransferData(anchorset = anchors, refdata = reference$cellTypes, dims = 1:30)
  obj <- Seurat::AddMetaData(obj, metadata = predictions)

  counts <- table(obj$predicted.id)
  obj <- subset(obj, predicted.id %in% names(counts[counts >= min_cells]))
  labels <- unique(obj$predicted.id)
  Seurat::Idents(obj) <- factor(obj$predicted.id,
                                levels = c(intersect(celltype_order, labels), sort(setdiff(labels, celltype_order))))
  obj
}

#' Whether the Lin et al. 2020 data (cached object or raw 10x files) is set up
lin2020_available <- function() {
  file.exists(analysis_path("Han2020_PDAmet.Rdata")) ||
    (nzchar(Sys.getenv("TOLDC_LIN2020_DIR")) && nzchar(Sys.getenv("TOLDC_LIN2020_META")))
}

#' Query datasets used in the paper
#'
#' `load` returns the raw query object given lazily loaded Lin et al. 2020
#' immune cells and the in-house annotated spleens; `available` says whether
#' its inputs are present.
external_datasets <- list(
  list(
    id = "Lin2020_primary", figure = "SupplFigS7", sample_col = "ID", available = lin2020_available,
    load = function(lin2020, spleen) {
      obj <- subset(lin2020(), Type %in% "Primary")
      # Remove mast cells
      subset(obj, cells = colnames(obj)[obj@assays$RNA$data["CPA3", ] <= 0.2])
    }
  ),
  list(
    id = "Chijimatsu2022_primary", figure = "SupplFigS7", sample_col = "Patient", scale_all = FALSE,
    available = function() file.exists(input_path("PDA_REF", "immunesubset.rds")),
    load = function(lin2020, spleen) {
      subset(readRDS(input_path("PDA_REF", "immunesubset.rds")), Type %in% "Tumor")
    }
  ),
  list(
    id = "Lin2020_liver_metastasis", figure = "SupplFigS8", sample_col = "ID", available = lin2020_available,
    resolution = seq(0, 2, 0.1),
    load = function(lin2020, spleen) subset(lin2020(), Type %in% "Metastasis" & Site == "Liver")
  ),
  list(
    id = "Normal_spleen", figure = "Fig1E_right", sample_col = "orig.ident", available = function() TRUE,
    resolution = seq(0, 2, 0.1),
    load = function(lin2020, spleen) subset(spleen(), sample %in% "Normal")
  )
)

#' Run one [external_datasets] entry and write its panels
run_external_dataset <- function(ds, reference, lin2020, spleen) {
  message("External dataset: ", ds$id)
  obj <- transfer_labels(
    ds$load(lin2020, spleen), reference,
    resolution = ds$resolution, scale_all = !identical(ds$scale_all, FALSE)
  )
  prefix <- paste0(ds$figure, "_", ds$id)
  save_table(as.data.frame(table(sample = obj[[ds$sample_col, drop = TRUE]], cell_type = obj$predicted.id)),
             paste0(prefix, "_cell_counts.tsv"))

  if (ds$figure != "Fig1E_right") {
    save_pdf(composition_bar_plot(obj, ds$sample_col, "predicted.id"),
             paste0(prefix, "_composition.pdf"), width = 6, height = 5)
    save_pdf(composition_bar_plot(obj, ds$sample_col, "predicted.id", y = "n", celltypes = mdc_subsets),
             paste0(prefix, "_mDC_counts.pdf"), width = 6, height = 5)
    for (spec in external_lr_specs) {
      save_pdf(lr_violin(obj, spec$genes, spec$cells, ylim = if (is.null(spec$ylim)) c(1, 5.5) else spec$ylim),
               paste0(prefix, "_", spec$name, "_violin.pdf"), width = 6, height = 6)
    }
  }

  present <- intersect(mdc_subsets, unique(obj$predicted.id))
  if (length(present) == 0) {
    message("  no mDC subclusters transferred; skipping signature panels")
    return(invisible(obj))
  }
  scored <- add_signature_score(subset(obj, predicted.id %in% mdc_subsets))
  if (ds$figure != "Fig1E_right") {
    save_pdf(feature_umap(scored, c("tol_DC", "IDO1"), pt_size = 0.5),
             paste0(prefix, "_signature_umap.pdf"), width = 10, height = 5)
  }
  save_pdf(signature_violin(scored, present), paste0(prefix, "_signature_violin.pdf"), width = 4, height = 5)
  invisible(obj)
}
