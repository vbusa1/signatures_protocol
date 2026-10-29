library(bluster)
library(cluster)
library(ggplot2)

set.seed(40)

# ==============================================================================
# STEP 2: Find Optimal Resolution
# ==============================================================================

# Compute nearest-neighbor graph
seu_merged <- FindNeighbors(seu_merged, dims = 1:15)

# Extract low-dimensional coordinates for distance calculations
pca_mat  <- Embeddings(seu_merged, "pca")[, 1:15]
dist_mat <- dist(pca_mat)

resolutions <- seq(0.1, 2, by = 0.1)
metrics     <- data.frame()

for (res in resolutions) {
  seu_merged <- FindClusters(seu_merged, 
                             resolution = res, 
                             verbose = FALSE)
  clusters   <- Idents(seu_merged)
  
  # Calculate quantitative metrics
  purity <- mean(neighborPurity(pca_mat, clusters)$purity)
  sil    <- mean(silhouette(as.numeric(clusters), dist_mat)[, 3])
  
  metrics <- rbind(metrics, data.frame(Resolution = res, 
                                       Purity = purity, 
                                       Silhouette = sil))
}

ggplot(metrics, aes(x = Resolution)) +
  geom_line(aes(y = Purity, color = "Purity")) +
  geom_line(aes(y = Silhouette, color = "Silhouette")) +
  theme_bw() + 
  labs(y = "Score", color = NULL)

ggsave("figs/02_metrics.pdf", width = 5, height = 3.5)

# Visualize cluster stability across evaluated resolutions
library(clustree)

clustree(
  seu_merged,
  prefix      = "RNA_snn_res.",
  node_colour = "grey80",
  edge_width  = 0.5,  
  show_axis   = TRUE
) + guides(edge_alpha = "none")

ggsave("figs/02_clustree.pdf", width = 7, height = 9)


# Compare optimal vs over-clustered resolutions side-by-side
seu_merged <- FindClusters(seu_merged, 
                           resolution = 0.5, 
                           verbose = FALSE)
seu_merged$cell_type_clust <- Idents(seu_merged)

# Contrast against over-clustered resolution
seu_merged <- FindClusters(seu_merged,
                           resolution = 2.0, 
                           verbose = FALSE)
seu_merged$res_overclustered <- Idents(seu_merged)

DimPlot(seu_merged, group.by = c("cell_type_clust", 
                                 "res_overclustered")) &
  theme_void()

ggsave("figs/02_dimplot.pdf", width = 8, height = 4)



# Check cross-species proportion within each cluster
ggplot(seu_merged@meta.data, 
       aes(x = cell_type_clust, 
           fill = species)) +
  geom_bar(position = "fill") +
  theme_classic() +
  labs(x = "Cluster", 
       y = "Proportion", 
       fill = NULL) +
  guides(fill = "none")

ggsave("figs/02_prop1.pdf", width = 3, height = 4)


ggplot(seu_merged@meta.data, 
       aes(x = res_overclustered, 
           fill = species)) +
  geom_bar(position = "fill") +
  theme_classic() +
  labs(x = "Cluster", 
       y = "Proportion", 
       fill = NULL)

ggsave("figs/02_prop2.pdf", width = 5, height = 4)





