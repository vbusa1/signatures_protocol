library(Seurat)
library(dplyr)
library(ggplot2)

set.seed(40)

# ==============================================================================
# STEP 1: Rename  -> Intersect -> Subset & Merge -> Visualize
# ==============================================================================

# Read in files
ortho    <- read.csv("ortholog_table.csv", 
                     check.names = FALSE, 
                     stringsAsFactors = FALSE)
seu_list <- readRDS("multispecies_seurat_list.rds")

# 1. Rename features to reference names and add species metadata
for (sp in names(seu_list)) {
  counts <- GetAssayData(seu_list[[sp]], layer = "counts")
  
  # using column 'reference_name' as gene name across all species
  map <- setNames(ortho$reference_name, ortho[[sp]])
  rownames(counts) <- unname(map[rownames(counts)])
  
  seu_list[[sp]] <- CreateSeuratObject(counts = counts, 
                                       meta.data = seu_list[[sp]]@meta.data)
  seu_list[[sp]]$species <- sp
}

# 2. Identify shared 1-to-1 orthologs across species and subset
shared_genes <- Reduce(intersect, lapply(seu_list, rownames))
seu_list <- lapply(seu_list, subset, features = shared_genes)

# 3. Merge into unified Seurat object
seu_merged <- merge(seu_list[[1]], 
                    y = seu_list[-1], 
                    add.cell.ids = names(seu_list))
seu_merged <- JoinLayers(seu_merged)


seu_merged <- seu_merged |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 2000) |>
  ScaleData() |>
  RunPCA() |>
  RunUMAP(dims = 1:15, 
          n.neighbors = 30, 
          min.dist = 0.3) |>
  FindNeighbors(dims = 1:15) |>
  FindClusters(resolution = 0.2)

DimPlot(seu_merged, group.by = c("seurat_clusters", 
                                 "cell_type", 
                                 "species")) &
  theme_void()

ggsave("figs/01_dimplots.pdf", width = 12, height = 4)

# DimPlot(seu_merged, reduction = "umap") # 3 clusters, corresponding broadly to cell type
# DimPlot(seu_merged, group.by = "cell_type")  # Shows 3 connected, distinct lobes
# DimPlot(seu_merged, group.by = "species")    # Shows all 4 species evenly intermixed inside each lobe

# ==============================================================================
# Optional Harmony correction
# ==============================================================================

library(harmony)

seu_merged <- RunHarmony(seu_merged, group.by.vars = "species")

# Re-run UMAP using the harmonized embeddings
seu_merged <- RunUMAP(seu_merged, 
                      reduction = "harmony", 
                      dims = 1:15)|>
  FindNeighbors(dims = 1:15) |>
  FindClusters(resolution = 0.2)

DimPlot(seu_merged, group.by = "species") &
  theme_void()



