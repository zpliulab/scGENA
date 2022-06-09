
# Step 6: Imputation 
# We used SAVER method

# Input data: differential genes matrix [ 2169g  1472 cells]
load("Differential_expression_Xin_data.RData")
# subset the expression for selected DE genes

DE_genes <- DE_genes[,1]
DE_mat <- xin_nonNorm[DE_genes,]
dim(DE_mat)  #2169 1472

# save original data for later comparison
library(data.table)
fwrite(x = DE_mat, file = "output/before_imputation.csv", row.names = TRUE)

library(SAVER)

DE_Mat_saver <- saver(DE_mat, 
                      ncores = 12, 
                      estimates.only = TRUE)
dim(DE_Mat_saver) #2169 1472
DE_saver <- as.data.frame(DE_Mat_saver)

save(DE_mat, DE_saver, file = "Differential Mat before & aft imp.RData")


#plot Pseudotime
# https://satijalab.org/signac/articles/monocle.html
# https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

# remotes::install_github('satijalab/seurat-wrappers')
# devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))
# load libraries
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)


#plot trajectories colored by pseudotime 
CellSets <- as.cell_data_set(xinse_meta)
CellSets <- cluster_cells(CellSets,
                          reduction_method = "UMAP")
CellSets <- learn_graph(CellSets, use_partition = TRUE)

CellSets <- order_cells(CellSets, 
                        reduction_method = "UMAP", root_cells = NULL)

tiff("figures/psedotime trejctories1.tiff", units="in", width=9, height=5, res=300)

plot_cells(
  cds = CellSets,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, label_roots = F)
dev.off()

tiff("figures/celltype trejctories 1.tiff", units="in", width=9, height=5, res=300)

plot_cells(
  cds = CellSets,
  color_cells_by = "cell_type1",
  show_trajectory_graph = TRUE, label_roots = F)
dev.off()