# Load necessary libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(hdf5r)

# Define the raw data path
data_path <- "/content/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5"

# Read the 10X Genomics data
nsclc.sparse.m <- Read10X_h5(data_path, 
                             use.names = TRUE, 
                             unique.features = TRUE)

# Extract counts matrix
cts <- nsclc.sparse.m[[1]]

# Create a Seurat object
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, 
                                       project = "NSCLC", 
                                       min.cells = 3, 
                                       min.features = 200)

# Calculate mitochondrial percentage
nsclc.seurat.obj$mito.percent <- PercentageFeatureSet(object = nsclc.seurat.obj, pattern = "^MT-")

# Filter cells based on quality control metrics
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mito.percent < 5)

# Normalize the data
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)

# Identify highly variable features
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Scale the data
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = rownames(nsclc.seurat.obj))

# Perform PCA
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(nsclc.seurat.obj))

# Doublet detection using DoubletFinder
sweep.res.list_nsclc <- paramSweep(nsclc.seurat.obj, PCs = 1:20, sct = FALSE)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)

# Extract the optimal pK value
pK <- bcmvn_nsclc %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) %>%
  as.numeric(as.character(.[[1]]))

# Estimate homotypic doublet proportion
annotations <- nsclc.seurat.obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)

# Adjust the expected number of doublets
nExp_poi <- round(0.076 * nrow(nsclc.seurat.obj@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
nsclc.seurat.obj <- doubletFinder(nsclc.seurat.obj, 
                                  PCs = 1:20, 
                                  pN = 0.25, 
                                  pK = pK, 
                                  nExp = nExp_poi.adj, 
                                  reuse.pANN = FALSE, 
                                  sct = FALSE)

# Filter out detected doublets
nsclc.seurat.obj.dblt.adj <- subset(nsclc.seurat.obj, subset = DF.classifications_0.25_0.08_1642 == "Singlet")

# Continue with downstream analysis using filtered Seurat object
# Find neighbors and clusters
nsclc.seurat.obj.dblt.adj <- FindNeighbors(nsclc.seurat.obj.dblt.adj, dims = 1:15)
nsclc.seurat.obj.dblt.adj <- FindClusters(nsclc.seurat.obj.dblt.adj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

# Run UMAP for dimensionality reduction
nsclc.seurat.obj.dblt.adj <- RunUMAP(nsclc.seurat.obj.dblt.adj, dims = 1:20)

# Save the final Seurat object
saveRDS(nsclc.seurat.obj.dblt.adj, file = "nsclc_seurat_obj_dblt_adj.rds")

# Print completion message
cat("Single-cell RNA-seq analysis pipeline completed successfully.")


