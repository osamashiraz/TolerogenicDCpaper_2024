# Download the data needed to regenerate the figures:
#   - the integrated spleen object, Azimuth reference, cached intermediates
#     and the modified CellChat database used for the paper (Zenodo)
#   - the Tabula Sapiens spleen dataset (cellxgene; Suppl. Fig. S2 and
#     rebuilding the Azimuth reference)
# Files already present with the expected checksum are skipped.
#
#   Rscript scripts/00_download_data.R

source("scripts/_setup.R")
options(timeout = 6 * 3600)

zenodo_record <- Sys.getenv("TOLDC_ZENODO_RECORD", "24099265")
zenodo_url <- function(file) sprintf("https://zenodo.org/records/%s/files/%s?download=1", zenodo_record, file)
ts_spleen_id <- "cbab13ec-d454-4d3b-8ab0-86715b140d83"
azimuth_dir <- c("References", "TS_spleen", "azimuth", "reference")

downloads <- list(
  list(url = zenodo_url("Integrated_bf_CCA_Res0.6_short.Rdata"),
       dest = input_path("Integrated_bf_CCA_Res0.6_short.Rdata"), md5 = "7cc0c2c25fbf09917a0a1e96a380c63a"),
  list(url = zenodo_url("TS_spleen_azimuth_ref.Rds"),
       dest = do.call(input_path, as.list(c(azimuth_dir, "ref.Rds"))), md5 = "926a0e0fa84e703741000ed0fc551c5c"),
  list(url = zenodo_url("TS_spleen_azimuth_idx.annoy"),
       dest = do.call(input_path, as.list(c(azimuth_dir, "idx.annoy"))), md5 = "3a20f3e1b92ed7ab90d41bd020d30cf8"),
  list(url = zenodo_url("SpleenAnnotated_clean.Rdata"),
       dest = analysis_path("SpleenAnnotated_clean.Rdata"), md5 = "ada06824b6afd34276c12dc6a003654f"),
  list(url = zenodo_url("PDA_spleen.Rdata"),
       dest = analysis_path("PDA_spleen.Rdata"), md5 = "91a6030e089cbfc76bfd0e0608280b78"),
  list(url = zenodo_url("PDA_spleen_DC_subset.Rdata"),
       dest = analysis_path("PDA_spleen_DC_subset.Rdata"), md5 = "59e113fb8c665b3a5548c0e97ae8e1fe"),
  list(url = zenodo_url("PDA_spleen_cellchat.Rdata"),
       dest = analysis_path("PDA_spleen_cellchat.Rdata"), md5 = "24fff6dcf2ef1505679cd83be7319a04"),
  list(url = zenodo_url("mDC_subclusters_vs_eachother.tsv"),
       dest = analysis_path("mDC_subclusters_vs_eachother.tsv"), md5 = "8bcf435922d5498048b3c4d68ba6118c"),
  list(url = zenodo_url("mDC_subclusters_vs_all_others_markers.tsv"),
       dest = analysis_path("mDC_subclusters_vs_all_others_markers.tsv"), md5 = "608840adbc3a4582adba7252c8e1c4f3"),
  list(url = zenodo_url("mDC_clusters_vs_all_others_markers.tsv"),
       dest = analysis_path("mDC_clusters_vs_all_others_markers.tsv"), md5 = "8a83b233a32e06149b32192ce34b86fa"),
  list(url = zenodo_url("CellChatDB.human_AAT.rda"),
       dest = "data/CellChatDB.human_AAT.rda", md5 = "484fe679bab4f9accf1906d10a75b9b4"),
  list(url = sprintf("https://datasets.cellxgene.cziscience.com/%s.rds", ts_spleen_id),
       dest = input_path("References", "TS_spleen", paste0(ts_spleen_id, ".rds")), md5 = NA)
)

md5_ok <- function(path, md5) file.exists(path) && (is.na(md5) || unname(tools::md5sum(path)) == md5)

for (d in downloads) {
  if (md5_ok(d$dest, d$md5)) {
    message("Already present: ", d$dest)
    next
  }
  message("Downloading ", d$dest)
  dir.create(dirname(d$dest), recursive = TRUE, showWarnings = FALSE)
  utils::download.file(d$url, d$dest, mode = "wb")
  if (!md5_ok(d$dest, d$md5)) stop("Checksum mismatch for ", d$dest)
}
