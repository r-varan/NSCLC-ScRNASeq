# NSCLC-ScRNASeq
Personal sc-RNASeq Project that uncovers Tumour Micro-Environment (TME) via cell-subtype clustering

# Dependencies
R
tidyverse
Seurat
ggplot2
DoubletFinder
hdf5r

# Input Data
The raw input data file 20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5 was derived from 10X Genomics website. It is not stored in a repository due to large file size. For more information, feel free to contact me.

# Output
The processed Seurat object with cell type clustering should be saved in the results/ directory as nsclc_seurat_obj_dblt_adj.rds.

# Detailed Steps

1. Data Loading
The raw scRNA-seq data is loaded from an HDF5 file using the Read10X_h5 function.

2. Quality Control
Cells are filtered based on the number of detected features and the percentage of mitochondrial genes to ensure high-quality data for analysis.

3. Normalization
The data is normalized using the NormalizeData function to make gene expression levels comparable across cells.

4. Feature Selection
Highly variable features are identified using the FindVariableFeatures function, which are used in downstream analyses.

5. Scaling
The data is scaled using the ScaleData function to standardize gene expression values.

6. Dimensionality Reduction
Principal Component Analysis (PCA) is performed using the RunPCA function to reduce the dimensionality of the data and highlight the most important features.

7. Clustering
Cells are clustered using a graph-based clustering approach with the FindNeighbors and FindClusters functions.

8. Visualization
The clusters are visualized using UMAP with the RunUMAP function, allowing for the identification of distinct cell populations.

9. Doublet Detection
Doublets are detected using the DoubletFinder package, and the data is adjusted to remove these artifacts.

10. Save Results
The final Seurat object, with doublets removed, is saved as an RDS file for future use.
