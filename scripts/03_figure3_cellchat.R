# Fig. 3A-F: CellChat communication between mregDCs (mDC_1) and other immune
# cells in the PDAC spleens. Also Suppl. Figs. S4-S6 and Suppl. Tables S6-S7.

source("scripts/_setup.R")

pda <- get_pdac()
pda <- merge_dc_subset_labels(pda, get_mdc(pda))
pda_cc <- prepare_cellchat_input(pda)
rm(pda)
cellchat <- get_cellchat(pda_cc)

# Suppl. Tables S6-S7: significant interactions -----------------------------
mdc1_interactions <- significant_interactions(cellchat, sources = "CD141+ mDC_1")
save_table(mdc1_interactions, "SupplTable_S6_cellchat_mDC1_significant_interactions.tsv")
save_table(significant_interactions(cellchat), "SupplTable_S7_cellchat_all_significant_interactions.tsv")

# Fig. 3A-B ----------------------------------------------------------------
save_pdf(mdc1_interaction_heatmap(cellchat), "Fig3A_mDC1_interaction_weights_counts.pdf", width = 9, height = 4,
         png = TRUE)
save_pdf(mdc1_top_interactions_heatmap(cellchat, mdc1_interactions),
         "Fig3B_mDC1_top_interactions_heatmap.pdf", width = 8, height = 9, png = TRUE)

# Fig. 3C-F and Suppl. Figs. S4-S5: chord diagram + expression per LR pair ---
for (spec in lr_specs) {
  save_pdf(function() draw_lr_pair(cellchat, pda_cc, spec),
           paste0(spec$figure, "_", spec$name, ".pdf"), width = 8, height = 7)
}

# Suppl. Fig. S6: T cell - macrophage interactions -------------------------
t_and_macrophages <- c("Macrophage", grep("T", levels(Idents(pda_cc)), value = TRUE))
pairs <- list(c("CD86", "CTLA4"), c("ICAM1", "ITGAL"), c("HLA-DMA", "CD4"), c("HLA-DMB", "CD4"),
              c("HLA-DRB1", "CD4"), c("IFNG", "IFNGR1"), c("IFNG", "IFNGR2"))
save_pdf(function() for (genes in pairs) print(lr_violin(pda_cc, genes, t_and_macrophages, ylim = c(1, 5))),
         "SupplFigS6_Tcell_macrophage_violins.pdf", width = 7, height = 6)
