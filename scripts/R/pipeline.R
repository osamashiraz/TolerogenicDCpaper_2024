#' Pipeline stages with on-disk caching
#'
#' Each heavy step is cached in `02-DataAnalysis/` under the file and object
#' names used for the paper. By default a stage loads its cache; with
#' `options(toldc.recompute = TRUE)` (the scripts' `--recompute` flag) every stage
#' is recomputed from its inputs and its cache overwritten.

#' Load a cached object, or compute and cache it
#'
#' @param path A `.tsv` table, or an `.Rdata` file holding one object called
#'   `name`.
#' @param compute Function producing the object.
cached <- function(path, name, compute, recompute = getOption("toldc.recompute", FALSE)) {
  is_tsv <- grepl("\\.tsv$", path)
  if (!recompute && file.exists(path)) {
    message("Loading ", path)
    if (is_tsv) return(utils::read.delim(path))
    env <- new.env()
    load(path, envir = env)
    return(env[[name]])
  }
  if (!recompute) message(path, " not found; computing it")
  result <- compute()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is_tsv) {
    utils::write.table(result, path, row.names = FALSE, quote = FALSE, sep = "\t")
  } else {
    assign(name, result)
    save(list = name, file = path)
  }
  result
}

#' Stage 1: Tabula Sapiens spleen Azimuth reference directory
get_azimuth_reference <- function(recompute = getOption("toldc.recompute", FALSE)) {
  ref_dir <- input_path("References", "TS_spleen", "azimuth", "reference")
  if (recompute || !file.exists(file.path(ref_dir, "ref.Rds"))) {
    raw <- readRDS(input_path("References", "TS_spleen", "cbab13ec-d454-4d3b-8ab0-86715b140d83.rds"))
    save_azimuth_reference(build_azimuth_reference(raw), ref_dir)
  }
  ref_dir
}

#' Stage 2: all 7 spleens, Azimuth-annotated, with `cellTypes` labels
get_spleen <- function(recompute = getOption("toldc.recompute", FALSE)) {
  spleen <- cached(analysis_path("SpleenAnnotated_clean.Rdata"), "SpleenObj_clean", recompute = recompute,
    compute = function() {
      env <- new.env()
      load(input_path("Integrated_bf_CCA_Res0.6_short.Rdata"), envir = env)
      obj <- env$joined_integ_bf_0.6_CCA_short
      obj$orig.ident[obj$orig.ident %in% "pdac_BulkFrozen"] <- "pdac_55"
      obj <- annotate_with_azimuth(obj, get_azimuth_reference())

      gene_map <- new.env()
      load(external_path("biomart_gene_map"), envir = gene_map)
      ensembl_to_symbol(obj, gene_map$geneMap)
    })
  label_cell_types(spleen)
}

#' Stage 3a: PDAC spleens (cell types without subcluster labels)
get_pdac <- function(recompute = getOption("toldc.recompute", FALSE)) {
  cached(analysis_path("PDA_spleen.Rdata"), "pda", recompute = recompute,
         compute = function() subset_pdac(get_spleen()))
}

#' Stage 3b: subclustered CD141+ mDCs, labeled mDC_1-3
get_mdc <- function(pda, recompute = getOption("toldc.recompute", FALSE)) {
  dc <- cached(analysis_path("PDA_spleen_DC_subset.Rdata"), "dc_SpleenObj", recompute = recompute,
               compute = function() subcluster_mdc(pda))
  assign_dc_subset_labels(dc)
}

#' Stage 4: CellChat object
get_cellchat <- function(pda_cc, recompute = getOption("toldc.recompute", FALSE)) {
  cached(analysis_path("PDA_spleen_cellchat.Rdata"), "cellchat", recompute = recompute,
         compute = function() build_cellchat(pda_cc))
}

#' Stage 5: immune cells of the Lin et al. 2020 dataset
#'
#' Immune clusters were identified at resolution 0.1 by marker expression
#' (PTPRC, CD3E, MS4A1, LYZ, CD74, ...). The cluster IDs belong to the cached
#' clustering; re-check them if `Han2020_PDAmet.Rdata` is recomputed.
get_lin2020_immune <- function(recompute = getOption("toldc.recompute", FALSE), immune_clusters = c(4, 9)) {
  obj <- cached(analysis_path("Han2020_PDAmet.Rdata"), "hanPDA_obj", recompute = recompute,
                compute = function() load_lin2020(external_path("lin2020_dir"), external_path("lin2020_meta")))
  subset(obj, RNA_snn_res.0.1 %in% immune_clusters)
}
