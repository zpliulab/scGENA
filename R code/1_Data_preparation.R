
# single cell Gene co-Expression Network Analysis (scGENA)

#in this pipeline we illustrate all the process for analysing single cell data
# uing R code
# we used Xin data to perform this pipeline   GEO ID: GSE81608
# set directory
setwd("F:/PhD/diabetes_datasets/scGena/scGENA_code")

# Input data


### step 1: Input the Data
xindat <- read.table("Data\GSE81608_human_islets_rpkm.txt", header = T)
dim(xindat)    # 39851  1601

genes <- read.csv("Data\\human_gene_annotation.csv", header = T)  
View(genes)

rownames(xindat) <- genes[,2]
xindat <- xindat[,2:ncol(xindat)]
dim(xindat)   # 39851  1600

dense.size <- object.size(as.matrix(xindat))
dense.size

sparse.size <- object.size(xindat)
sparse.size
dense.size/sparse.size

# 1.1: data exploration 

counts_per_cell <- Matrix::colSums(xindat)
cat("counts per cell: ", counts_per_cell[1:5], "\n") ## counts for first 5 cells

counts_per_gene <- Matrix::rowSums(xindat)
cat("counts per gene: ", counts_per_gene[1:5], "\n")  ## counts for first 5 genes

genes_per_cell <- Matrix::colSums(xindat > 0) # count gene only if it has non-zero reads mapped.
cat("counts for non-zero genes: ", genes_per_cell[1:5])  ## counts for first 5 genes

#### cells_per_gene <- Matrix::?(counts>?) # only count cells where the gene is expressed
cells_per_gene <- Matrix::rowSums(xindat > 0) # only count cells where the gene is expressed
cat("count of cells with expressed genes: ", cells_per_gene)
dim(cells_per_gene)
dim(xindat)


hist(log10(counts_per_cell+1),main='counts per cell',col='dark grey')
dev.off()

hist(log10(genes_per_cell+1), main='genes per cell', col='dark grey')

plot(counts_per_cell, genes_per_cell, log='xy', col='dark grey')
title('counts vs genes per cell')

hist(log10(counts_per_gene+1), main='counts per genes', col='dark grey')

plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')


# 1.2: imput metadata

ann <- read.table("Data\\human_islet_cell_identity.txt", 
                  header = T, sep = "\t", stringsAsFactors = F)
dim(ann)   # 1600   14

rownames(ann) <- ann[,1]
rownames(ann) <- gsub(" ", "_", rownames(ann))
ann <- ann[,9:ncol(ann)]
colnames(ann)[length(colnames(ann))] <- "cell_type1"

# format cell type names
ann$cell_type1[ann$cell_type1 == "PP"] <- "gamma"
ann$cell_type1[ann$cell_type1 == "PP.contaminated"] <- "gamma.contaminated"

# # # #remove contaminated cells from cells
# ann<-ann[!(ann$cell_type1=="beta.contaminated" | ann$cell_type1=="alpha.contaminated" | ann$cell_type1=="delta.contaminated"| ann$cell_type1=="gamma.contaminated"),]
# 

#extract the celltypes
celltypes <- ann[,6]
celltypes <- matrix(celltypes)
colnames(celltypes) <- "cell_type1"

#combine celltypes to remove contaminated cells
xindat <- as.data.frame(t(xindat))
dim(xindat)  # 1600 39851

xindat <- cbind(xindat, celltypes)
dim(xindat)  # 1600 39852

# #remove contaminated cells
xindat<-xindat[!(xindat$cell_type1=="beta.contaminated" | xindat$cell_type1=="alpha.contaminated" | xindat$cell_type1=="delta.contaminated"| xindat$cell_type1=="gamma.contaminated"),]
dim(xindat)   #  1492 39852

xindat <- as.data.frame(t(xindat))
dim(xindat) #  39852  1492 
xindat <- xindat[-c(39852),]  #  39851  1492 

save(xindat, file = "Xin_data_for_preprocessing.RData")

# Next step prepocessing 