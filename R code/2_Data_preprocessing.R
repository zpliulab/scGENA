
load("Xin_data_for_preprocessing.RData")

# step 2: Preprocessing of the data

# load required libraries

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(ggthemes)
  library(dplyr)
  library(ggplot2)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(knitr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(data.table)
  library(RColorBrewer)
  library(MAST)
  library(cowplot)
  library(patchwork)
  library(limma)
  library(Biobase)
  library(marray)
  library(convert)
  library(Matrix)
})

#we use seurat for preprocessing analysis

## Initialize the Seurat object with the raw (non-normalized data).
xinseurat <- CreateSeuratObject(counts = xindat, project = "islet", 
                                min.cells = 3, 
                                min.features = 500)
dim(xinseurat)    # 29043  1492
levels(xinseurat)    # sample



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
xinseurat[["percent.mt"]] <- PercentageFeatureSet(xinseurat, pattern = "^MT")
dim(xinseurat) #29043  1492

# Visualize QC metrics as a violin plot
# pdf(file = "figures/QC_vilon_plot.pdf", wi = 10, he = 6)
# tiff("figures/QC_vilon.tiff", units="in", width=7, height=4, res=300)

VlnPlot(xinseurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# dev.off()


# Show QC metrics for the first 5 cells
head(xinseurat@meta.data, 5)


# pdf(file = "figures/QC feature scatters.pdf", wi = 12, he = 6)
# tiff("figures/QC features sactters.tiff", units="in", width=12, height=5, res=300)

plot1 <- FeatureScatter(xinseurat, 
                        feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(xinseurat, 
                        feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# dev.off()

#
xinseurat <- subset(xinseurat, 
                    subset = nFeature_RNA > 1500 & nFeature_RNA < 10000 & percent.mt < 1.5)

tiff("figures/QC features.tiff", units="in", width=12, height=5, res=300)

p1 <- FeatureScatter(xinseurat, 
                     feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(xinseurat, 
                     feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(p1, p2))
dev.off()



#data normalization
#normalized data are saved in [RNA]@data

xinseurat <- NormalizeData(xinseurat)

# xinseurat <- NormalizeData(xinseurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
xinseurat <- FindVariableFeatures(xinseurat, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(xinseurat), 10)
top10

# plot variable features with and without labels
# pdf(file = "figures/top10 variable features.pdf", wi = 12, he = 6)
tiff("figures/variable features.tiff", units="in", width=12, height=6, res=300)

plot1 <- VariableFeaturePlot(xinseurat)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2
dev.off()


## data scaling 
# The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(xinseurat)
xinseurat <- ScaleData(xinseurat, features = all.genes)
head(xinseurat[["RNA"]]@scale.data,5)
dim(xinseurat)  #29510  1579
levels(xinseurat)


## Perform dimensional reduction
xinseurat <- RunPCA(xinseurat, features = VariableFeatures(object = xinseurat))
# Examine and visualize PCA results a few different ways
print(xinseurat[["pca"]], dims = 1:5, nfeatures = 5)

# pdf(file = "figures/PCA linear dimensional reduction.pdf", wi = 12, he = 8)
# tiff("figures/Dim reduction.tiff", units="in", width=12, height=8, res=300)

VizDimLoadings(xinseurat, dims = 1:2, reduction = "pca")

# dev.off()

# tiff("figures/dimplot.tiff", units="in", width=12, height=5, res=300)

DimPlot(xinseurat, reduction = "pca")
# dev.off()

# tiff("figures/dimheatmap.tiff", units="in", width=10, height=6, res=300)

DimHeatmap(xinseurat, dims = 1, cells = 500, balanced = TRUE)
# dev.off()

# tiff("figures/dimheatmap 9.tiff", units="in", width=10, height=6, res=300)

DimHeatmap(xinseurat, dims = 1:9, cells = 500, balanced = TRUE)
# dev.off()





## Cluster the cells
xinseurat <- FindNeighbors(xinseurat, dims = 1:10)
xinseurat <- FindClusters(xinseurat, resolution = 0.7)  # resolution to incr/dec No. of clusters
## Look at cluster IDs of the first 5 cells
head(Idents(xinseurat), 5)
levels(xinseurat)  #



#UMAP
# xinseurat <- RunUMAP(xinseurat, dims = 1:10, verbose = F)
# tiff("figures/clusters UMAP.tiff", units="in", width=9, height=5, res=300)

xinseurat <- RunUMAP(xinseurat, dims = 1:10)

DimPlot(xinseurat, reduction = "umap")
# dev.off()


#tsne
# tiff("figures/tSNE plot clusters.tiff", units="in", width=9, height=5, res=300)

xinseurat <- RunTSNE(xinseurat, dims = 1:10)

DimPlot(xinseurat, reduction = "tsne")
DimPlot(xinseurat, reduction = "tsne", label = TRUE)
# dev.off()
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

save(xinseurat, file = "seurat_object_of_Xin_Data.RData")

# Next step 5: DE genes
