#' Tabula Sapiens spleen Azimuth reference (Methods: "Cell type identification")
#'
#' Follows the Azimuth reference-building helper script
#' (satijalab/azimuth-references, `export.R`).

#' Build an Azimuth reference from the Tabula Sapiens spleen object
#'
#' @param ref_obj Seurat object downloaded from cellxgene, with
#'   `cell_ontology_class`, `tissue` and `compartment` metadata.
#' @param downsample_n Cells kept per `cell_ontology_class`.
#' @return The `Azimuth::AzimuthReference()` object.
build_azimuth_reference <- function(ref_obj, downsample_n = 2000) {
  Seurat::Idents(ref_obj) <- "cell_ontology_class"
  ref_obj$cell_ontology_class <- as.character(ref_obj$cell_ontology_class)

  # The published reference was downsampled without a fixed seed; the seed
  # here makes rebuilds repeatable but not identical to it.
  set.seed(123)
  ref_obj <- subset(ref_obj, downsample = downsample_n)

  options(future.globals.maxSize = 1591289600)
  ref_obj <- Seurat::SCTransform(ref_obj)
  ref_obj <- Seurat::RunPCA(ref_obj)
  ref_obj <- Seurat::FindNeighbors(ref_obj, dims = 1:30)
  ref_obj <- Seurat::FindClusters(ref_obj)
  ref_obj <- Seurat::RunUMAP(ref_obj, dims = 1:30, return.model = TRUE)

  Azimuth::AzimuthReference(
    object = ref_obj,
    refUMAP = "umap",
    refDR = "pca",
    refAssay = "SCT",
    metadata = c("tissue", "cell_ontology_class", "compartment"),
    dims = 1:50,
    k.param = 31,
    reference.version = "1.0.0"
  )
}

#' Write an Azimuth reference in the directory layout `RunAzimuth()` expects
#'
#' @param azimuth_ref Output of [build_azimuth_reference()].
#' @param ref_dir Directory to write `ref.Rds` and `idx.annoy` into.
save_azimuth_reference <- function(azimuth_ref, ref_dir) {
  dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
  SeuratObject::SaveAnnoyIndex(
    object = azimuth_ref[["refdr.annoy.neighbors"]],
    file = file.path(ref_dir, "idx.annoy")
  )
  saveRDS(azimuth_ref, file = file.path(ref_dir, "ref.Rds"))
  invisible(ref_dir)
}
