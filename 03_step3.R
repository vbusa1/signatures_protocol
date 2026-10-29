library(Seurat)

set.seed(40)

# ==============================================================================
# STEP 3: Identify Signature Genes
# ==============================================================================

# =======
# Part A
# =======

seu_list <- SplitObject(seu_merged, split.by = "species")

markers_list <- list()

for (sp in names(seu_list)) {
  # Set active identity to the evaluated cell types
  Idents(seu_list[[sp]]) <- "cell_type"
  
  # Calculate positive markers relative to species background
  markers_list[[sp]] <- FindAllMarkers(
    seu_list[[sp]], 
    only.pos = TRUE,
    logfc.threshold = 0.5,
    min.pct = 0.25,
    min.diff.pct = 0.10
  )
}


cell_types <- unique(seu_merged$cell_type)
conserved_markers <- list()

for (ct in cell_types) {
  genes_per_sp <- lapply(markers_list, function(m) {
    subset(m, cluster == ct)$gene
  })
  
  # Retain only genes identified as markers in ALL evaluated species
  conserved_markers[[ct]] <- Reduce(intersect, genes_per_sp)
}

# ## non-conserved markers:
# CT_3 :
#   [1] "GENE-1128" "GENE-1168" "GENE-1184" "GENE-1292" "GENE-1434" "GENE-1450" "GENE-1485" "GENE-1520" "GENE-1605" "GENE-1624" "GENE-171" 
# [12] "GENE-172"  "GENE-173"  "GENE-174"  "GENE-1740" "GENE-175"  "GENE-176"  "GENE-177"  "GENE-178"  "GENE-180"  "GENE-181"  "GENE-182" 
# [23] "GENE-1825" "GENE-183"  "GENE-1838" "GENE-184"  "GENE-186"  "GENE-187"  "GENE-188"  "GENE-189"  "GENE-190"  "GENE-1905" "GENE-191" 
# [34] "GENE-193"  "GENE-194"  "GENE-195"  "GENE-196"  "GENE-197"  "GENE-1971" "GENE-198"  "GENE-199"  "GENE-200"  "GENE-201"  "GENE-202" 
# [45] "GENE-203"  "GENE-204"  "GENE-205"  "GENE-206"  "GENE-208"  "GENE-209"  "GENE-210"  "GENE-223"  "GENE-292"  "GENE-299"  "GENE-801" 
# [56] "GENE-804"  "GENE-807"  "GENE-811"  "GENE-818"  "GENE-824"  "GENE-828"  "GENE-843"  "GENE-854"  "GENE-868"  "GENE-870"  "GENE-874" 
# [67] "GENE-875"  "GENE-884"  "GENE-887"  "GENE-890"  "GENE-906"  "GENE-912"  "GENE-915"  "GENE-922"  "GENE-933"  "GENE-939"  "GENE-954" 
# [78] "GENE-964"  "GENE-968"  "GENE-970"  "GENE-974"  "GENE-977"  "GENE-981"  "GENE-989"  "GENE-993"  "GENE-995"  "GENE-999" 
# 
# CT_1 :
#   [1] "GENE-100"  "GENE-1007" "GENE-101"  "GENE-102"  "GENE-1028" "GENE-103"  "GENE-104"  "GENE-105"  "GENE-106"  "GENE-107"  "GENE-108" 
# [12] "GENE-109"  "GENE-110"  "GENE-111"  "GENE-112"  "GENE-113"  "GENE-114"  "GENE-1140" "GENE-1142" "GENE-1143" "GENE-1146" "GENE-115" 
# [23] "GENE-116"  "GENE-1179" "GENE-118"  "GENE-119"  "GENE-120"  "GENE-121"  "GENE-1217" "GENE-123"  "GENE-125"  "GENE-126"  "GENE-127" 
# [34] "GENE-128"  "GENE-1281" "GENE-129"  "GENE-130"  "GENE-1384" "GENE-1470" "GENE-1530" "GENE-1535" "GENE-1556" "GENE-1659" "GENE-1662"
# [45] "GENE-1666" "GENE-168"  "GENE-1715" "GENE-176"  "GENE-18"   "GENE-1854" "GENE-1869" "GENE-1873" "GENE-1879" "GENE-1884" "GENE-1986"
# [56] "GENE-210"  "GENE-216"  "GENE-393"  "GENE-398"  "GENE-408"  "GENE-422"  "GENE-428"  "GENE-432"  "GENE-434"  "GENE-443"  "GENE-447" 
# [67] "GENE-448"  "GENE-449"  "GENE-456"  "GENE-459"  "GENE-465"  "GENE-469"  "GENE-471"  "GENE-474"  "GENE-476"  "GENE-487"  "GENE-490" 
# [78] "GENE-497"  "GENE-513"  "GENE-519"  "GENE-526"  "GENE-537"  "GENE-550"  "GENE-553"  "GENE-555"  "GENE-558"  "GENE-562"  "GENE-565" 
# [89] "GENE-569"  "GENE-573"  "GENE-576"  "GENE-577"  "GENE-579"  "GENE-592"  "GENE-91"   "GENE-92"   "GENE-93"   "GENE-94"   "GENE-95"  
# [100] "GENE-96"   "GENE-97"   "GENE-98"   "GENE-99"  
# 
# CT_2 :
#   [1] "GENE-1054" "GENE-1096" "GENE-1119" "GENE-125"  "GENE-131"  "GENE-132"  "GENE-133"  "GENE-134"  "GENE-135"  "GENE-136"  "GENE-137" 
# [12] "GENE-1375" "GENE-138"  "GENE-139"  "GENE-140"  "GENE-141"  "GENE-142"  "GENE-143"  "GENE-144"  "GENE-145"  "GENE-146"  "GENE-147" 
# [23] "GENE-148"  "GENE-149"  "GENE-1491" "GENE-150"  "GENE-151"  "GENE-152"  "GENE-153"  "GENE-154"  "GENE-155"  "GENE-1551" "GENE-156" 
# [34] "GENE-157"  "GENE-158"  "GENE-159"  "GENE-1591" "GENE-160"  "GENE-1606" "GENE-161"  "GENE-162"  "GENE-163"  "GENE-164"  "GENE-1641"
# [45] "GENE-165"  "GENE-166"  "GENE-1667" "GENE-167"  "GENE-168"  "GENE-169"  "GENE-170"  "GENE-32"   "GENE-400"  "GENE-601"  "GENE-609" 
# [56] "GENE-624"  "GENE-637"  "GENE-639"  "GENE-642"  "GENE-643"  "GENE-658"  "GENE-660"  "GENE-678"  "GENE-694"  "GENE-700"  "GENE-707" 
# [67] "GENE-708"  "GENE-714"  "GENE-721"  "GENE-722"  "GENE-726"  "GENE-728"  "GENE-735"  "GENE-745"  "GENE-753"  "GENE-758"  "GENE-775" 
# [78] "GENE-791"  "GENE-799"

FeaturePlot(
  seu_merged,
  features = c(
    "GENE-194", 
    "GENE-97",
    "GENE-161"  
    ),
  pt.size = 1.5,
  split.by = "species",
  reduction = "umap"
) & theme_void() & NoLegend()
ggsave("figs/03_dimplot_nonconsmarkers.pdf", width = 8, height = 8) 



## marker genes that are DEGs
# "GENE-982" "GENE-445" "GENE-618" "GENE-671" "GENE-732" "GENE-796"

FeaturePlot(
  seu_merged,
  features = c("GENE-982", "GENE-618", "GENE-671"),
  pt.size = 1.5,
  split.by = "species",
  reduction = "umap"
) & theme_void() & NoLegend()
ggsave("figs/03_dimplot_DEGmarkers.pdf", width = 8, height = 8) 


# =======
# Part B
# =======


library(DESeq2)

species_list <- unique(seu_merged$species)
cell_types <- unique(seu_merged$cell_type)
comb_list <- combn(species_list, 2, 
                   simplify = FALSE) # All unique species pairs

tdsq_list <- list()

for (ct in cell_types) {
  # 1. Subset target cell type
  seu_ct <- subset(seu_merged, cell_type == ct)
  
  # 2. Define grouping variable for pseudo-bulk aggregation
  # paste the species and replicate metadata columns via "."
    seu_ct$pb_group <- paste(seu_ct$species, 
                             seu_ct$replicate, sep = ".")
  
  # 3. Aggregate raw counts to pseudo-bulk matrix
  pb_counts <- AggregateExpression(seu_ct, 
                                   group.by = "pb_group", 
                                   slot = "counts")$RNA
  
  # 4. Construct sample metadata (colData)
  # re-split by "." to re-claim the species
  col_species <- sapply(strsplit(colnames(pb_counts), 
                                 ".", fixed = TRUE), `[`, 1)
  col_data <- data.frame(
    species = factor(col_species),
    row.names = colnames(pb_counts)
  )
  
  # 5. Build DESeq2 object and fit model
  # incorporate additional covariates such as batch if needed
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_counts),
    colData = col_data,
    design = ~ species
  )
  
  tdsq_list[[ct]] <- DESeq(dds, quiet = TRUE)
}


non_deg_genes <- list()

for (ct in cell_types) {
  tdsq <- tdsq_list[[ct]]
  all_genes <- rownames(tdsq)
  stable_genes <- all_genes
  
  # Evaluate stability across all pairwise species comparisons
  for (comb in comb_list) {
    res <- DESeq2::results(
      tdsq, 
      lfcThreshold = 1.5,
      altHypothesis = "lessAbs",
      contrast = c("species", comb[1], comb[2]),
      alpha = 0.05
    )
    
    # Retain genes that are significantly STABLE (padj < 0.05)
    # Ensure we drop NAs which can occur for lowly expressed genes
    res_sig <- subset(res, padj < 0.05 & !is.na(padj))
    non_deg <- rownames(res_sig)
    
    stable_genes <- intersect(stable_genes, non_deg)
  }
  
  non_deg_genes[[ct]] <- stable_genes
}



# =======
# Part C
# =======

conserved_signatures <- list()

for (ct in cell_types) {
  markers_A <- conserved_markers[[ct]]
  non_degs_B <- non_deg_genes[[ct]]
  
  candidate_signature <- intersect(markers_A, non_degs_B)
  
  if (length(candidate_signature) > 0) {
    conserved_signatures[[ct]] <- candidate_signature
  }
}



FeaturePlot(
  seu_merged,
  features = c(
    "GENE-72" , 
    "GENE-598", 
    "GENE-42" 
    ),
  pt.size = 1.5,
  split.by = "species",
  reduction = "umap"
)
ggsave("figs/03_dimplot_sig.pdf", width = 8, height = 8) 




library(ComplexHeatmap)

seu_merged$sp_ct <- paste(seu_merged$species, seu_merged$cell_type, sep = ".")
avg_exp <- AverageExpression(seu_merged, group.by = "sp_ct", slot = "data")$RNA



# Optional
all_sig_genes <- unique(unlist(conserved_signatures))



mat <- t(avg_exp[all_sig_genes, ])

row_info     <- strsplit(rownames(mat), "\\.")
species_vec  <- sapply(row_info, `[`, 1)
celltype_vec <- sapply(row_info, `[`, 2)

row_order    <- order(celltype_vec, species_vec)
mat          <- mat[row_order, ]
celltype_vec <- celltype_vec[row_order]
species_vec  <- species_vec[row_order]

mat_scaled <- scale(mat)
rownames(mat_scaled) <- species_vec

# Optional
mat_scaled[mat_scaled > 2.5]  <- 2.5
mat_scaled[mat_scaled < -2.5] <- -2.5

heatmap_composite <- Heatmap(
  mat_scaled,
  name = "Z-score",
  row_split = celltype_vec,
  cluster_rows = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  show_column_dend = FALSE
)

pdf("figs/03_heatmap.pdf", width = 7, height = 4)
draw(heatmap_composite)
dev.off()


