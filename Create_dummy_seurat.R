library(Seurat)
library(Matrix)

set.seed(42)

# --- 1. Dataset Parameters ---
species_list <- c("Species_A", "Species_B", "Species_C", "Species_D")
cell_types   <- c("CT_1", "CT_2", "CT_3")
n_cells_per_species <- 250 # 1,000 cells total
total_cells <- length(species_list) * n_cells_per_species
n_genes <- 2000

# --- 2. Setup Cell Metadata & Trajectory Coordinate (t) ---
# t_val represents cell progression from central origin (t = 0) to lobe tip (t = 1)
t_val <- rbeta(total_cells, shape1 = 2, shape2 = 2) 

meta_data <- data.frame(
  cell_id   = paste0("Cell_", 1:total_cells),
  species   = rep(species_list, each = n_cells_per_species),
  cell_type = rep(sample(cell_types, total_cells, replace = TRUE))
)
rownames(meta_data) <- meta_data$cell_id

# --- 3. Construct Baseline Expression Matrix ---
# Base expression vector (log-normal distribution)
base_gene_means <- rgamma(n_genes, shape = 2, scale = 0.3)
expected_counts <- matrix(rep(base_gene_means, total_cells), nrow = n_genes, ncol = total_cells)

# Add baseline noise
expected_counts <- expected_counts * matrix(rlnorm(n_genes * total_cells, meanlog = 0, sdlog = 0.2), nrow = n_genes)

# --- 4. Inject Pinwheel Trajectory Signals ---

# A. Central Hub Anchor Genes (Genes 211 - 400)
# Expressed highest at t = 0 across ALL cell types & species to bind the pinwheel center together
for (i in 1:total_cells) {
  expected_counts[211:400, i] <- expected_counts[211:400, i] * (1.5 + 2.5 * (1.0 - t_val[i]))
}

# B. Lobe Continuity Drivers (Genes 401 - 1000)
# Smoothly increase along t for each specific cell type (drives continuous lobe formation)
for (i in 1:total_cells) {
  ct <- meta_data$cell_type[i]
  t  <- t_val[i]
  if (ct == "CT_1") expected_counts[401:600, i]  <- expected_counts[401:600, i]  * (1.0 + 3.0 * t)
  if (ct == "CT_2") expected_counts[601:800, i]  <- expected_counts[601:800, i]  * (1.0 + 3.0 * t)
  if (ct == "CT_3") expected_counts[801:1000, i] <- expected_counts[801:1000, i] * (1.0 + 3.0 * t)
}

# C. Conserved Signature Genes (30 per CT across ALL species; Genes 1 - 90)
# Exponential gain toward the tips (t = 1) ensuring high FC for FindAllMarkers intersection
for (i in 1:total_cells) {
  ct <- meta_data$cell_type[i]
  t  <- t_val[i]
  if (ct == "CT_1") expected_counts[1:30, i]   <- expected_counts[1:30, i]   * (1.0 + 6.0 * (t^1.5))
  if (ct == "CT_2") expected_counts[31:60, i]  <- expected_counts[31:60, i]  * (1.0 + 6.0 * (t^1.5))
  if (ct == "CT_3") expected_counts[61:90, i]  <- expected_counts[61:90, i]  * (1.0 + 6.0 * (t^1.5))
}

# D. Species-Specific Markers (10 per CT per Species; Genes 91 - 210)
# Moderate expression gain within each species' subgroup without fragmenting the lobe
idx <- 91
for (ct in cell_types) {
  for (sp in species_list) {
    genes_sp <- idx:(idx + 9)
    for (i in 1:total_cells) {
      if (meta_data$cell_type[i] == ct && meta_data$species[i] == sp) {
        t <- t_val[i]
        expected_counts[genes_sp, i] <- expected_counts[genes_sp, i] * (1.0 + 3.5 * t)
      }
    }
    idx <- idx + 10
  }
}

# --- 5. Generate Final Sparse Counts & Process Seurat Object ---
# Add sequencing depth variance
size_factors <- rlnorm(total_cells, meanlog = 0, sdlog = 0.2)
expected_counts <- sweep(expected_counts, 2, size_factors, "*")

# Negative binomial sampling to simulate real scRNA-seq overdispersion
count_matrix <- matrix(rnbinom(length(expected_counts), size = 5, mu = expected_counts), nrow = n_genes)
rownames(count_matrix) <- paste0("GENE_", 1:n_genes)
colnames(count_matrix) <- meta_data$cell_id

count_matrix <- as(count_matrix, "sparseMatrix")

# Build and compute dimensional reduction
seu <- CreateSeuratObject(counts = count_matrix, meta.data = meta_data)
seu$pseudotime <- t_val
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, 
                            selection.method = "vst", 
                            nfeatures = 2000, 
                            verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:15, n.neighbors = 30, min.dist = 0.3, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu, dims = 1:15, resolution = .2)
DimPlot(seu, reduction = "umap", label = TRUE) # 3 clusters, corresponding broadly to cell type
DimPlot(seu, group.by = "cell_type")  # Shows 3 connected, distinct lobes
DimPlot(seu, group.by = "species")    # Shows all 4 species evenly intermixed inside each lobe

# Save the object locally
saveRDS(seu, "multispecies_seurat.rds")


