library(Seurat)
library(Matrix)

set.seed(40)

# ==============================================================================
# PART 1: Generate Master Multi-Species Seurat Object
# ==============================================================================

# --- 1. Dataset Parameters ---
species_list <- c("Species-A", "Species-B", "Species-C", "Species-D")
cell_types   <- c("CT_1", "CT_2", "CT_3")
n_cells_per_species <- 250 # 1,000 total cells across 4 species
total_cells <- length(species_list) * n_cells_per_species
n_genes <- 2000

# --- 2. Create Cell Metadata & Developmental Trajectory (t) ---
# t_val represents continuous cell progression from origin (t = 0) to tip (t = 1)
t_val <- rbeta(total_cells, shape1 = 2, shape2 = 2) 

meta_data <- data.frame(
  cell_id   = paste0("Cell_", 1:total_cells),
  species   = rep(species_list, each = n_cells_per_species),
  cell_type = sample(cell_types, total_cells, replace = TRUE)
)
rownames(meta_data) <- meta_data$cell_id

# --- 3. Construct Baseline Expression Matrix ---
# Draw baseline mean expression levels from a Gamma distribution
base_gene_means <- rgamma(n_genes, shape = 2, scale = 0.3)
expected_counts <- matrix(rep(base_gene_means, total_cells), nrow = n_genes, ncol = total_cells)

# Multiply by log-normal noise to simulate cell-to-cell biological variance
expected_counts <- expected_counts * matrix(rlnorm(n_genes * total_cells, meanlog = 0, sdlog = 0.2), nrow = n_genes)

# --- 4. Inject Trajectory & Cell-Type Expression Signals ---

# A. Central Hub Anchor Genes (Genes 211 - 400): High expression at origin (t = 0) across all species
for (i in 1:total_cells) {
  expected_counts[211:400, i] <- expected_counts[211:400, i] * (1.5 + 2.5 * (1.0 - t_val[i]))
}

# B. Lobe Continuity Drivers (Genes 401 - 1000): Increase along t for specific cell types
for (i in 1:total_cells) {
  ct <- meta_data$cell_type[i]
  t  <- t_val[i]
  if (ct == "CT_1") expected_counts[401:600, i]  <- expected_counts[401:600, i]  * (1.0 + 3.0 * t)
  if (ct == "CT_2") expected_counts[601:800, i]  <- expected_counts[601:800, i]  * (1.0 + 3.0 * t)
  if (ct == "CT_3") expected_counts[801:1000, i] <- expected_counts[801:1000, i] * (1.0 + 3.0 * t)
}

# C. Conserved Markers (Genes 1 - 90): Strong signal at trajectory tips (t = 1) across all species
for (i in 1:total_cells) {
  ct <- meta_data$cell_type[i]
  t  <- t_val[i]
  if (ct == "CT_1") expected_counts[1:30, i]   <- expected_counts[1:30, i]   * (1.0 + 6.0 * (t^1.5))
  if (ct == "CT_2") expected_counts[31:60, i]  <- expected_counts[31:60, i]  * (1.0 + 6.0 * (t^1.5))
  if (ct == "CT_3") expected_counts[61:90, i]  <- expected_counts[61:90, i]  * (1.0 + 6.0 * (t^1.5))
}

# D. Species-Specific Markers (Genes 91 - 210): Expression gains unique to specific species/cell type pairs
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

# --- 5. Generate Sparse Count Matrix & Build Base Seurat Object ---
# Add sequencing depth variance across cells
size_factors <- rlnorm(total_cells, meanlog = 0, sdlog = 0.2)
expected_counts <- sweep(expected_counts, 2, size_factors, "*")

# Sample counts using Negative Binomial distribution to simulate scRNA-seq overdispersion
count_matrix <- matrix(rnbinom(length(expected_counts), size = 5, mu = expected_counts), nrow = n_genes)
rownames(count_matrix) <- paste0("GENE-", 1:n_genes)
colnames(count_matrix) <- meta_data$cell_id

# Convert to memory-efficient sparse matrix
count_matrix <- as(count_matrix, "sparseMatrix")



# Assign 2 biological replicates (Rep1, Rep2) per species and cell type
# (needed to get DESeq2 to run)
meta_data$replicate <- "Rep1"
for (sp in species_list) {
  for (ct in cell_types) {
    idx <- which(meta_data$species == sp & meta_data$cell_type == ct)
    if (length(idx) > 0) {
      # Alternate Rep1 and Rep2 deterministically
      meta_data$replicate[idx] <- rep(c("Rep1", "Rep2"), length.out = length(idx))
    }
  }
}

seu <- CreateSeuratObject(counts = count_matrix, meta.data = meta_data)
saveRDS(seu, "multispecies_seurat.rds")




# ==============================================================================
# PART 2: Disassemble into Species List, Rename Features & Generate Ortholog Table
# ==============================================================================

species_list <- unique(seu$species)

# 1. Construct 1-to-1 Conserved Ortholog Mapping (GENE-1 to GENE-2000)
conserved_orthos <- data.frame(
  ortholog_id    = paste0("ORTHO-", 1:2000),
  reference_name = paste0("GENE-", 1:2000),
  stringsAsFactors = FALSE
)

for (sp in species_list) {
  sp_letter <- gsub("Species-", "", sp)
  conserved_orthos[[sp]] <- paste0("GENE-", 1:2000, "-", sp_letter)
}

# 2. Split Master Object into Species-Specific Seurat Objects & Add Non-Conserved Genes
seu_list <- list()
unshared_rows <- list()

for (sp in species_list) {
  sp_letter <- gsub("Species-", "", sp)
  
  # Extract cells belonging to current species
  sp_cells <- rownames(seu@meta.data[seu@meta.data$species == sp, ])
  seu_sp   <- subset(seu, cells = sp_cells)
  
  # Remove species column from metadata so downstream pipeline can re-add it cleanly
  seu_sp$species <- NULL
  
  # Extract raw count matrix and append species suffix to feature names
  counts_sp <- GetAssayData(seu_sp, layer = "counts")
  rownames(counts_sp) <- paste0(rownames(counts_sp), "-", sp_letter)
  
  # Generate 200 species-unique (non-conserved) genes
  n_extra <- 200
  extra_counts <- matrix(
    rnbinom(n_extra * ncol(counts_sp), size = 5, mu = 2), 
    nrow = n_extra, 
    ncol = ncol(counts_sp)
  )
  extra_gene_names <- paste0("GENE-", 1:n_extra, sp_letter, "1")
  rownames(extra_counts) <- extra_gene_names
  colnames(extra_counts) <- colnames(counts_sp)
  
  # Combine conserved and species-unique matrices
  combined_counts <- rbind(counts_sp, as(extra_counts, "sparseMatrix"))
  
  # Create individual species Seurat object
  seu_list[[sp]] <- CreateSeuratObject(
    counts = combined_counts, 
    meta.data = seu_sp@meta.data
  )
  
  # Record non-conserved rows for ortholog table (check.names = FALSE preserves hyphens)
  unshared_sp <- data.frame(
    ortholog_id    = paste0("UNSHARED-", sp_letter, "-", 1:n_extra),
    reference_name = paste0("GENE-", 1:n_extra, sp_letter),
    check.names    = FALSE,
    stringsAsFactors = FALSE
  )
  
  for (s_col in species_list) {
    unshared_sp[[s_col]] <- if (s_col == sp) extra_gene_names else NA_character_
  }
  
  unshared_rows[[sp]] <- unshared_sp
}

# Combine conserved and non-conserved ortholog mappings into master table
ortholog_table <- rbind(conserved_orthos, do.call(rbind, unshared_rows))

# 3. Save Final Outputs
write.csv(ortholog_table, "ortholog_table.csv", row.names = FALSE)
saveRDS(seu_list, "multispecies_seurat_list.rds")
