#' CellChat ligand-receptor analysis of the PDAC spleens (Fig. 3)

#' Ligand-receptor pairs shown as chord diagrams and violin plots
#'
#' `signaling`/`pair` identify the pair in CellChat; `cells` are the cell
#' types in the violin panel, or `"interacting"` to use every cell type that
#' receives the signal from mDC_1. `figure` is the output file prefix.
lr_specs <- list(
  list(figure = "Fig3C", name = "IDO1-AHR", signaling = "IDO1", pair = "IDO1-AHR",
       genes = c("IDO1", "AHR"),
       cells = c(mdc_subsets, "T_reg cell", "Stem cell", "Monocyte", "Macrophage",
                 "Innate Lymphoid cell", "CD4+ T_naive cell")),
  list(figure = "Fig3D", name = "CCL19-CCR7", signaling = "CCL", pair = "CCL19_CCR7",
       genes = c("CCL19", "CCR7"),
       cells = c(mdc_subsets, "B_mem cell", "B_naive cell", "Innate Lymphoid cell",
                 "CD4+ T_naive cell", "CD4+ T_mem cell", "T_reg cell")),
  list(figure = "Fig3E", name = "NECTIN2-TIGIT", signaling = "NECTIN", pair = "NECTIN2_TIGIT",
       genes = c("NECTIN2", "TIGIT"), ylim = c(1, 5),
       cells = c(mdc_subsets, "T_reg cell", "CD4+ T_mem cell", "CD4+ T_naive cell", "CD8+ T cell",
                 "CD8+ T_mem cell", "CD8+ T_naive cell", "NK cell")),
  list(figure = "Fig3F", name = "CLEC2D-KLRB1", signaling = "CLEC", pair = "CLEC2D_KLRB1",
       genes = c("CLEC2D", "KLRB1"),
       cells = c(mdc_subsets, "B_mem cell", "B_naive cell", "Innate Lymphoid cell",
                 "CD4+ T_naive cell", "CD4+ T_mem cell", "T_reg cell")),
  list(figure = "SupplFigS4", name = "AHR-IDO1", signaling = "AHR", pair = "AHR-IDO1",
       genes = c("IDO1", "AHR"),
       cells = c(mdc_subsets, "T_reg cell", "Stem cell", "Monocyte", "Macrophage",
                 "Innate Lymphoid cell", "CD4+ T_naive cell", "CD4+ T_mem cell")),
  list(figure = "SupplFigS5", name = "CD86-CTLA4", signaling = "CD86", pair = "CD86_CTLA4",
       genes = c("CD86", "CTLA4"),
       cells = c(mdc_subsets, "T_reg cell", "Monocyte", "Macrophage", "B_mem cell", "CD4+ T_mem cell")),
  list(figure = "SupplFigS5", name = "CDH1-KLRG1", signaling = "CDH1", pair = "CDH1_KLRG1",
       genes = c("CDH1", "KLRG1"), cells = "interacting"),
  list(figure = "SupplFigS5", name = "PD-L1-PD-1", signaling = "PD-L1", pair = "CD274_PDCD1",
       genes = c("CD274", "PDCD1"), cells = "interacting")
)

#' Prepare the PDAC object for CellChat
#'
#' Drops cell types with fewer than 15 cells (erythrocytes, platelets,
#' neutrophils) and orders identities for plotting.
#'
#' @param pda PDAC object with subcluster labels in `cellTypes`.
prepare_cellchat_input <- function(pda) {
  pda <- subset(pda, cellTypes %notin% c("Erythrocyte", "Platelet", "Neutrophil"))
  pda$samples <- pda$orig.ident
  levels <- intersect(celltype_order, unique(pda$cellTypes))
  pda$cellTypes <- factor(pda$cellTypes, levels = levels)
  Seurat::Idents(pda) <- pda$cellTypes
  Seurat::DefaultAssay(pda) <- "RNA"
  pda
}

#' Run CellChat inference with the modified ligand-receptor database
#'
#' @param pda_cc Output of [prepare_cellchat_input()].
#' @param db_path `.rda` holding the modified CellChatDB as `db.AAT`.
build_cellchat <- function(pda_cc, db_path = "data/CellChatDB.human_AAT.rda") {
  cellchat <- CellChat::createCellChat(object = pda_cc, group.by = "cellTypes")
  db <- new.env()
  load(db_path, envir = db)
  cellchat@DB <- db$db.AAT

  cellchat <- CellChat::subsetData(cellchat)
  future::plan("multisession", workers = 4)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(
    cellchat, type = "truncatedMean", trim = 0.1,
    distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
    contact.dependent = TRUE, contact.range = 100
  )
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  CellChat::aggregateNet(cellchat)
}

#' Significant (p <= 0.05) ligand-receptor interactions as a table
#'
#' @param sources Restrict to signals sent by these cell types (`NULL` = all).
significant_interactions <- function(cellchat, sources = NULL) {
  CellChat::netVisual_bubble(cellchat, sources.use = sources, return.data = TRUE, thresh = 0.05)$communication
}

#' Fig. 3A: interaction weights and counts from mDC_1 to every cell type
mdc1_interaction_heatmap <- function(cellchat, from = "CD141+ mDC_1") {
  weights <- t(as.matrix(cellchat@net$weight[from, ]))
  counts <- t(as.matrix(cellchat@net$count[from, ]))
  rownames(weights) <- "Interaction\nWeights"
  rownames(counts) <- "Interaction\nCounts"
  ramp <- function(m) {
    circlize::colorRamp2(c(0, stats::quantile(m)[3:5]), c("#2C7BB6", "#ABD9E9", "#FDAE61", "#D7191C"))
  }

  ComplexHeatmap::Heatmap(
    weights, name = "Interaction Weights", col = ramp(weights), height = grid::unit(1, "cm"),
    heatmap_legend_param = list(direction = "horizontal"), border = TRUE
  ) %v% ComplexHeatmap::Heatmap(
    counts, name = "Interaction Counts", col = ramp(counts), height = grid::unit(1, "cm"),
    heatmap_legend_param = list(direction = "horizontal"), border = TRUE,
    column_names_gp = grid::gpar(fontsize = 12, col = celltype_cols[colnames(weights)])
  )
}

#' Fig. 3B: strongest interactions from mDC_1 to its top interactors
#'
#' Top interactors are cell types receiving a total interaction weight above
#' `min_weight` from mDC_1; interactions shown are those with a scaled
#' communication score (the `prob` column of `netVisual_bubble()`) of at
#' least `min_score` in one of them. Rows are alphabetical.
#'
#' @param mdc1_interactions [significant_interactions()] for mDC_1.
mdc1_top_interactions_heatmap <- function(cellchat, mdc1_interactions,
                                          min_weight = 0.25, min_score = 0.25) {
  df <- mdc1_interactions
  df$label <- as.character(df$interaction_name_2)
  df$label <- gsub("-", " => ", df$label)
  df$label <- gsub("HLA => ", "HLA-", df$label)
  df$label <- gsub("  =>  ", " => ", df$label)

  score <- reshape2::dcast(df, label ~ target, value.var = "prob", fill = 0)
  rownames(score) <- score$label
  score <- as.matrix(score[, -1])

  weights <- cellchat@net$weight["CD141+ mDC_1", ]
  top_interactors <- names(weights[weights > min_weight])
  mat <- score[, intersect(colnames(score), top_interactors)]
  mat <- mat[rowSums(mat >= min_score) > 0, ]

  cell_group <- ifelse(grepl("B_", colnames(mat)), "B cell",
                ifelse(grepl("CD141+", colnames(mat), fixed = TRUE), "mDC",
                ifelse(grepl("NK", colnames(mat)), "NK cell",
                ifelse(grepl("T_|T cell", colnames(mat)), "T cell", "Other"))))
  cell_group <- factor(cell_group, levels = c("mDC", "T cell", "NK cell", "B cell", "Other"))

  ComplexHeatmap::Heatmap(
    mat, name = "Interaction Scores", border = TRUE,
    col = circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
                               c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = ifelse(grepl("HLA", rownames(mat)), "Antigen\nPresentation", "Other"),
    column_split = cell_group,
    column_names_gp = grid::gpar(fontsize = 12, col = celltype_cols[colnames(mat)])
  )
}

#' Chord diagram, CellChat gene-expression violins and per-cell-type violins
#' for one [lr_specs] entry, drawn on the current device
draw_lr_pair <- function(cellchat, pda_cc, spec) {
  cells <- spec$cells
  if (identical(cells, "interacting")) {
    prob <- cellchat@net$prob[, , spec$pair]["CD141+ mDC_1", ]
    cells <- c("CD141+ mDC_1", names(prob[prob > 0]))
  }
  color_use <- celltype_cols[colnames(cellchat@net$weight)]

  CellChat::netVisual_individual(
    cellchat, signaling = spec$signaling, pairLR.use = spec$pair, layout = "chord",
    cell.order = levels(Seurat::Idents(pda_cc)), color.use = color_use
  )
  print(CellChat::plotGeneExpression(
    cellchat, signaling = spec$signaling, enriched.only = TRUE, type = "violin", color.use = color_use
  ))
  print(lr_violin(pda_cc, spec$genes, cells, ylim = if (is.null(spec$ylim)) c(1, 5.5) else spec$ylim))
}
