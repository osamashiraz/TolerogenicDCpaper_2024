# mDC subclusters and key LR pairs in other datasets, mapped onto the
# in-house PDAC spleens by label transfer:
#   Fig. 1E (right)   in-house normal spleens
#   Suppl. Fig. S7    primary PDAC: Lin et al. 2020, Chijimatsu et al. 2022
#   Suppl. Fig. S8    PDAC liver metastases: Lin et al. 2020
#
# Run a subset with e.g. `--dataset=Normal_spleen,Chijimatsu2022_primary`.

source("scripts/_setup.R")

pda <- get_pdac()
pda <- merge_dc_subset_labels(pda, get_mdc(pda))

# Load the larger shared inputs once, and only if a selected dataset needs them.
lin2020 <- local({
  obj <- NULL
  function() {
    if (is.null(obj)) obj <<- get_lin2020_immune()
    obj
  }
})
spleen <- local({
  obj <- NULL
  function() {
    if (is.null(obj)) obj <<- get_spleen()
    obj
  }
})

selected <- sub("^--dataset=", "", grep("^--dataset=", commandArgs(trailingOnly = TRUE), value = TRUE))
datasets <- external_datasets
if (length(selected) > 0) {
  datasets <- Filter(function(ds) ds$id %in% strsplit(selected, ",")[[1]], datasets)
}

failed <- character()
for (ds in datasets) {
  if (!ds$available()) {
    message("Skipping ", ds$id, ": input data not found (see README, step 4)")
    next
  }
  tryCatch(
    run_external_dataset(ds, pda, lin2020, spleen),
    error = function(e) {
      message("  FAILED ", ds$id, ": ", conditionMessage(e))
      failed <<- c(failed, ds$id)
    }
  )
}
if (length(failed) > 0) stop("Datasets not completed: ", paste(failed, collapse = ", "))
