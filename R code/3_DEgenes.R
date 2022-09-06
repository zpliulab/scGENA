
## step 5: Differential Gene expression (DGE) analysis


# load seurat object of thedata

load("seurat_object_of_Xin_Data.RData")

#run required data

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
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


#we use find all markers by seurat
xinseurat.markers <- FindAllMarkers(xinseurat, only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.25)
xinseurat.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
dim(xinseurat.markers)    #2169    7

#identify top10 differentiated genes
top10 <- xinseurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Visualise several genes 

VlnPlot(xinseurat, features = c("IAPP", "ID2"))

VlnPlot(xinseurat, features = c("TTR", "INS", "TM4SF4", "PPY", "GCG", "RBP4", "IAPP", "FABP5","SPINK1"))

FeaturePlot(xinseurat, features = c("TTR", "INS", "TM4SF4", "PPY", "GCG", "RBP4", "IAPP", "SST", "SPINK1"))

RidgePlot(xinseurat, features = c("PPY", "ID2"))

DotPlot(xinseurat, features = c("PPY", "ID2", "IAPP", "GCG"))

DoHeatmap(xinseurat, features = top10$gene, label = TRUE)#+ NoLegend()

DoHeatmap(subset(xinseurat, downsample = 10), features = features, size = 3)


# ## find Cluster that match the cell types of data
# resolution can be changed to fit the number of cluster
xinseurat <- FindNeighbors(xinseurat, dims = 1:10)
xinseurat <- FindClusters(xinseurat, resolution = 0.1)  # resolution to incr/dec No. of clusters
head(Idents(xinseurat), 5)

# # # #remove contaminated cells from cells
# ann <- ann[!(ann$cell_type1=="beta.contaminated" | ann$cell_type1=="alpha.contaminated" | ann$cell_type1=="delta.contaminated"| ann$cell_type1=="gamma.contaminated"),]

# Identify cell type ----------------------------------------------------------------
celltype <- ann$cell_type1
celltype[!duplicated(celltype)]
which(celltype == "alpha")    # 886
#which(celltype == "alpha.contaminated")    # 60
which(celltype == "beta")    # 472
#which(celltype == "beta.contaminated")    # 31
which(celltype == "delta")    # 49
#which(celltype == "delta.contaminated")    # 9
which(celltype == "gamma")    # 85   PP
#which(celltype == "gamma.contaminated")    # 8

#check the cells number
levels(xinseurat)

#new.cluster.ids <- c("alpha", "beta", "gamma", "alpha.contaminated", "delta", "beta.contaminated", "delta.contaminated", "gamma.contaminated")
new.cluster.ids <- c("alpha", "beta", "gamma", "delta")
names(new.cluster.ids) <- levels(xinseurat)
xinseurat <- RenameIdents(xinseurat, new.cluster.ids)
DimPlot(xinseurat, reduction = "umap", label = TRUE, pt.size = 0.1) #+ NoLegend()
DimPlot(xinseurat, reduction = "tsne", label = TRUE, pt.size = 0.1) #+ NoLegend()


DimPlot(xinseurat, reduction = "umap", label = TRUE)

DimPlot(xinseurat, reduction = "tsne", label = TRUE)

VlnPlot(xinseurat, features = "percent.mt", 
        # split.pl = "data",
)
VlnPlot(xinseurat, features = "nFeature_RNA", 
        # split.pl = "data",
)


# # Differential expression using MAST
# BiocManager::install("MAST")
library(MAST)

xinseura.MAST <-FindAllMarkers(xinseurat,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25,
                               test.use = "MAST") 
dim(xinseura.MAST)   



#combine metadata to seurat to extract genes in each cell-type
library(SeuratObject)
xinse_meta <- AddMetaData(xinseurat, ann)

#DE by cell types
xinse_meta.markers1 <- FindAllMarkers(xinse_meta, only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25)

#gen_alph <- DE_cell_markers[1:133,]    #133
#gen_beta <- DE_cell_markers[134:1403,] #1270
#gen_gam <- DE_cell_markers[1404:1748,]  #345
#gen_deta <- DE_cell_markers[1749:1919,] #171

#create matrix to genes
#gen_alph <- xin_nonNorm[gen_alph,]
#gen_beta <- xin_nonNorm[gen_beta,]
#gen_gam <- xin_nonNorm[gen_gam,]
#gen_deta <- xin_nonNorm[gen_deta,]



save(xinseurat, xinseurat.markers, xinse_meta, 
     gen_alph,gen_beta, gen_deta, gen_gam, 
     file = "Differential_expression_Xin_data.RData")

# Next step: Imputation
