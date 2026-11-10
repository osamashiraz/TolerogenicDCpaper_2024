# Fig. 1A-E: cell types of normal and PDAC spleens, CD141+ mDC subclusters,
# their markers and the tolerogenic signature. Also Suppl. Fig. S3 and
# Suppl. Tables S2-S5. (The normal-spleen panel of Fig. 1E is made by
# 04_external_validation.R.)

source("scripts/_setup.R")

spleen <- get_spleen()
pda <- get_pdac()
dc <- get_mdc(pda)
pda <- merge_dc_subset_labels(pda, dc)
spleen <- merge_dc_subset_labels(spleen, dc)

save_table(as.data.frame(table(sample = spleen$orig.ident, cell_type = spleen$cellTypes)),
           "cell_counts_per_sample.tsv")

# Fig. 1A ------------------------------------------------------------------
save_pdf(umap_plot(spleen, "cellTypes", celltype_cols, reduction = "ref.umap"),
         "Fig1A_umap_cell_types.pdf", width = 9, height = 6, png = TRUE)
save_pdf(umap_plot(spleen, "sample", condition_cols, reduction = "ref.umap",
                   split.by = "sample", pt_size = 0.4, alpha = 0.4),
         "Fig1A_umap_by_condition.pdf", width = 10, height = 5)
save_pdf(umap_plot(dc, "dc_subsets", celltype_cols, pt_size = 1, alpha = 1) + NoLegend() + ggtitle(" "),
         "Fig1A_inset_mDC_subclusters_umap.pdf", width = 5, height = 5)

# Fig. 1B ------------------------------------------------------------------
save_pdf(composition_bar_plot(pda, "orig.ident", "cellTypes"),
         "Fig1B_composition_PDAC.pdf", width = 5, height = 5)
normal <- subset(spleen, sample %in% "Normal")
normal_counts <- table(normal$cellTypes)
normal <- subset(normal, cellTypes %in% names(normal_counts[normal_counts >= 10]))
save_pdf(composition_bar_plot(normal, "orig.ident", "cellTypes"),
         "Fig1B_composition_Normal.pdf", width = 4, height = 5)

# Suppl. Fig. S3: mDC subclusters per PDAC spleen ---------------------------
save_pdf(composition_bar_plot(dc, "orig.ident", "dc_subsets", y = "n", celltypes = mdc_subsets),
         "SupplFigS3_mDC_subcluster_counts.pdf", width = 6, height = 5)

# Fig. 1C ------------------------------------------------------------------
save_pdf(feature_umap(dc, c("IDO1", "IDO2")), "Fig1C_IDO1_IDO2_umap.pdf", width = 10, height = 5, png = TRUE)

# Suppl. Tables S2-S4: DEG comparisons (a-c in the Methods) -----------------
# The published tables are loaded from 02-DataAnalysis/ unless --recompute;
# Wilcoxon results can shift slightly across Seurat versions.
each_other <- cached(analysis_path("mDC_subclusters_vs_eachother.tsv"),
                     compute = function() markers_subclusters_vs_each_other(dc))
vs_all <- cached(analysis_path("mDC_subclusters_vs_all_others_markers.tsv"),
                 compute = function() markers_subclusters_vs_all(pda))
whole <- cached(analysis_path("mDC_clusters_vs_all_others_markers.tsv"),
                compute = function() markers_whole_cluster_vs_rest(pda, sig_filter(each_other)))
save_table(each_other, "SupplTable_S2_mDC_subclusters_vs_each_other.tsv")
save_table(vs_all, "SupplTable_S3_mDC_subclusters_vs_all_other_cells.tsv")
save_table(whole, "SupplTable_S4_mDC_cluster_vs_all_other_cells.tsv")

# Fig. 1D / Suppl. Table S5 ------------------------------------------------
each_other_sig <- sig_filter(each_other)
fig1d_markers <- build_figure1d_markers(
  each_other_sig = each_other_sig,
  vs_all_sig = sig_filter(vs_all),
  whole_sig = subset(sig_filter(whole), gene %notin% each_other_sig$gene)
)
save_table(fig1d_markers, "SupplTable_S5_mDC_subcluster_top_markers.tsv")

expr <- FetchData(dc, vars = unique(fig1d_markers$gene), layer = "data")
set.seed(123)
save_pdf(
  Heatmap(
    t(expr), name = "Log Expression", col = rev(RColorBrewer::brewer.pal(n = 4, name = "RdYlBu")),
    row_split = fig1d_markers$cluster[match(colnames(expr), fig1d_markers$gene)],
    column_split = dc$dc_subsets, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
    show_column_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
    row_names_side = "right", row_names_gp = gpar(fontface = "bold"),
    row_gap = unit(1, "mm"), column_gap = unit(1, "mm"), row_title = " ", column_title = " ",
    border = TRUE, heatmap_legend_param = list(direction = "horizontal")
  ),
  "Fig1D_mDC_marker_heatmap.pdf", width = 7, height = 7, png = TRUE
)

# Fig. 1E: tolerogenic signature in the PDAC mDC subclusters ----------------
dc <- add_signature_score(dc)
save_pdf(feature_umap(dc, "tol_DC", limits = c(-1, 1)), "Fig1E_signature_umap.pdf", width = 5, height = 5)
save_pdf(signature_violin(dc), "Fig1E_signature_violin_PDAC.pdf", width = 4, height = 5, png = TRUE)
