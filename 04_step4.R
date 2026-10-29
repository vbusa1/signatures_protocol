library(AUCell)

# AUCell operates on raw ranks, making it immune to normalization differences
expr_matrix <- GetAssayData(seu_merged, layer = "counts")
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)

cells_AUC <- AUCell_calcAUC(conserved_signatures, cells_rankings,
                            aucMaxRank = ceiling(0.05 * nrow(expr_matrix)))
auc_mat <- getAUC(cells_AUC)

# Append these scores to Seurat metadata for plotting
auc_features <- c()
for (ct in names(conserved_signatures)) {
  feat_name <- paste0("AUC_", ct)
  seu_merged[[feat_name]] <- auc_mat[ct, ]
  auc_features <- c(auc_features, feat_name)
}


# Spatial mapping of signature weights (UMAP)
FeaturePlot(seu_merged, 
            features = auc_features)

VlnPlot(seu_merged, 
        features = auc_features, 
        group.by = "cell_type",
        pt.size = 0)
