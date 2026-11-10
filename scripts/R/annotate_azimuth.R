#' Reference mapping of the integrated spleen dataset with Azimuth

#' Annotate the integrated spleen object against the Tabula Sapiens reference
#'
#' Runs `Azimuth::RunAzimuth()`, then drops unassigned cells and ontology
#' classes with fewer than `min_cells` cells.
#'
#' @param spleen_obj Integrated (CCA) Seurat object of all 7 spleens.
#' @param reference_dir Directory written by [save_azimuth_reference()].
#' @param min_cells Minimum cells for an ontology class to be kept.
#' @return The annotated, filtered Seurat object.
annotate_with_azimuth <- function(spleen_obj, reference_dir, min_cells = 10) {
  options(future.globals.maxSize = 1891289600)
  set.seed(123)
  spleen_obj <- Azimuth::RunAzimuth(spleen_obj, assay = "RNA", reference = reference_dir)

  counts <- table(spleen_obj$predicted.cell_ontology_class)
  rare <- names(counts[counts < min_cells])
  message("Dropping rare classes (< ", min_cells, " cells): ", paste(rare, collapse = ", "))
  keep <- !is.na(spleen_obj$predicted.cell_ontology_class) &
    spleen_obj$predicted.cell_ontology_class %notin% rare
  subset(spleen_obj, cells = colnames(spleen_obj)[keep])
}

#' Rename Ensembl-ID features to gene symbols and re-process
#'
#' Features without a symbol are dropped. Seurat keys every layer by feature
#' name, so normalization, scaling, PCA, clustering and UMAP are redone after
#' the rename. The Azimuth reference UMAP (`ref.umap`) is kept.
#'
#' @param seurat_obj Seurat v5 object with Ensembl-ID features.
#' @param gene_map data.frame with Ensembl-ID row names and an
#'   `external_gene_name` column (biomaRt export).
#' @return The renamed, re-processed Seurat object.
ensembl_to_symbol <- function(seurat_obj, gene_map) {
  ens <- rownames(seurat_obj@assays$RNA@features@.Data)
  symbol <- gene_map[ens, "external_gene_name"]
  keep <- ens[!is.na(symbol)]

  seurat_obj <- Seurat::DietSeurat(
    seurat_obj, assays = "RNA", features = keep, dimreducs = c("ref.umap", "umap")
  )
  new_names <- gene_map[rownames(seurat_obj@assays$RNA@features@.Data), "external_gene_name"]
  rownames(seurat_obj@assays$RNA@features@.Data) <- new_names
  for (layer in intersect(c("counts", "data", "scale.data"), names(seurat_obj@assays$RNA@layers))) {
    rownames(seurat_obj@assays$RNA@layers[[layer]]) <- new_names
  }

  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- Seurat::RunPCA(seurat_obj)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:15)
  seurat_obj <- Seurat::FindClusters(seurat_obj)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:15)
  Seurat::UpdateSeuratObject(seurat_obj)
}

#' Attach the paper's short cell-type labels
#'
#' @param seurat_obj Object with a `predicted.cell_ontology_class` column.
#' @return The object with `cellTypes` set from [ontology_to_label].
label_cell_types <- function(seurat_obj) {
  seurat_obj$cellTypes <- unname(ontology_to_label[seurat_obj$predicted.cell_ontology_class])
  seurat_obj
}
