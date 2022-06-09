
# Step 5: Gene co-expression network analysis
# we build a networks based on pseudo-bulk cell types

# packages
library(Seurat)
library(Matrix)
library(tidyverse)
library(rliger)
library(dplyr)

# library(scWGCNA)

# celltyps="alpha", "beta", "gamma", "delta"

load("RData/seurat_obj_with_metadata_celltypes.RData")
load("RData/Betaseurat_betaMat_betaSaver.RData")
load("RData/beta preprocessin & imp.RData")

seurat_beta <- subset(xinse_meta, cell_type1 %in% c("beta"))
dim(seurat_beta) #29043   460

# save(xinse_meta, file = "RData/seurat_obj_with_metadata_celltypes.RData")


seurat_mark <- FindAllMarkers(seurat_beta, only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

genes_mark <- seurat_mark[,7]
beta_mat <- as.data.frame(as.matrix(seurat_beta@assays[["RNA"]]@data))
sub_mar <- beta_mat[genes_mark,]

#Imputation method, skip if not necessary
library(SAVER)
Beta_saver <- saver(sub_mar, 
                    ncores = 12, 
                    estimates.only = TRUE)

save(Beta_saver, seurat_beta, seurat_mark, file = "RData/beta preprocessin & imp.RData")

#run GCNs

library(tidyverse)
library(WGCNA)
library(flashClust)
enableWGCNAThreads()


datExpr <- as.data.frame(t(Beta_saver))
# 
# # only keep good genes:
datExpr <- datExpr[,goodGenes(datExpr)]
save(datExpr, datTraits, file = "RData/Beta_datExpre.RData")

load("RData/Beta_datExpre.RData")


# # Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,20, by=2));

# # Call the network topology analysis function for each set in turn
powerTable = list(
  data = pickSoftThreshold(
    datExpr,
    powerVector=powers,
    verbose = 100,
    networkType="signed",
    corFnc="bicor"
  )[[2]]
)

# Plot the results:
pdf("figures/Beta_softPower.pdf", height=10, width=18)

colors = c("blue", "red","black")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "mean connectivity",
             "Max connectivity");

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (col in 1:length(plotCols)){
  ylim[1, col] = min(ylim[1, col], powerTable$data[, plotCols[col]], na.rm = TRUE);
  ylim[2, col] = max(ylim[2, col], powerTable$data[, plotCols[col]], na.rm = TRUE);
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;

for (col in 1:length(plotCols)){
  plot(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
       xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
       main = colNames[col]);
  addGrid();

  if (col==1){
    text(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
         labels=powers,cex=cex1,col=colors[1]);
  } else
    text(powerTable$data[,1], powerTable$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[1]);
  if (col==1){
    legend("bottomright", legend = 'Betacells', col = colors, pch = 20) ;
  } else
    legend("topright", legend = 'Betacells', col = colors, pch = 20) ;
}
dev.off()

# Based on the soft power threshold that you have selected, now we can build the co-expression network. Consult the WGCNA documentation if you need help selecting a soft power value. Furthermore, you should read the function description for blockwiseConsensusModules carefully to select the different parameters, however the ones that I have chosen here have generally given good results on a variety of datasets.



softPower=8

nSets = 1
setLabels = 'beta'
shortLabels = setLabels

multiExpr <- list()
multiExpr[['beta']] <- list(data=datExpr)

checkSets(multiExpr) # check data size


# construct network
net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                              maxBlockSize = 1000, ## This should be set to a smaller size if the user has limited RAM
                              randomSeed = 12345,
                              corType = "bicor", # "bicor"bidweight mid correlation
                              power = softPower,
                              consensusQuantile = 0.3,
                              networkType = "signed",
                              TOMType = "unsigned",
                              TOMDenom = "min",
                              scaleTOMs = TRUE, scaleQuantile = 0.8,
                              sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                              useDiskCache = TRUE, chunkSize = NULL,
                              deepSplit = 4,
                              pamStage=FALSE,
                              detectCutHeight = 0.995, minModuleSize = 50,
                              mergeCutHeight = 0.2,
                              saveConsensusTOMs = TRUE,
                              consensusTOMFilePattern = "ConsensusTOM-block.%b.rda")



consMEs = net$multiMEs;
moduleLabels = net$colors;

# Convert the numeric labels to color labels
moduleColors = as.character(moduleLabels)
consTree = net$dendrograms[[1]];

# module eigengenes
MEs=moduleEigengenes(multiExpr[[1]]$data, colors = moduleColors, nPC=1)$eigengenes
MEs=orderMEs(MEs)
meInfo<-data.frame(rownames(datExpr), MEs)
colnames(meInfo)[1]= "SampleID"
# View(moduleColors)
# intramodular connectivity
KMEs<-signedKME(datExpr, MEs,outputColumnName = "kME",corFnc = "bicor")

# compile into a module metadata table
geneInfo=as.data.frame(cbind(colnames(datExpr),moduleColors, KMEs))

# how many modules did we get?
nmodules <- length(unique(moduleColors))

# merged gene symbol column
colnames(geneInfo)[1]= "GeneSymbol"
colnames(geneInfo)[2]= "Initially.Assigned.Module.Color"

# save info
write.csv(geneInfo,file=paste0('output/geneInfoSigned.csv'))
# write.csv(TOM.matrix,file=paste0('output/TOM_matrix.csv'))
# TOM_1 <- as.data.frame(consTomDS)
# write.csv(as.data.frame(consTomDS),file=paste0('output/Tom_cons.csv'))

PCvalues=MEs


# pdf("figures/SignedDendro.pdf",height=5, width=8)
plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = paste0("Beta lineage gene dendrogram and module colors"))
# dev.off()


library(igraph)
library(RColorBrewer)

load("ConsensusTOM-block.1.rda")

TOM.matrix = as.matrix(consTomDS);
uniquemodcolors = unique(geneInfo$Initially.Assigned.Module.Color);
# uniquemodcolors=uniquemodcolors[uniquemodcolors!= 'blue']


# pdf(paste0('figures/Beta ModuleNetworks_4.pdf'),height=9,width=10)
# tiff("figures/dimplot.tiff", units="in", width=12, height=9, res=300)

for (mod in uniquemodcolors)  {
  numgenesingraph = 25;
  numconnections2keep = 500;
  cat('module:',mod,'\n');
  geneInfo=geneInfo[geneInfo$GeneSymbol!="NA",]
  colind = which(colnames(geneInfo)==paste('kME',mod,sep=''));
  rowind = which(geneInfo$Initially.Assigned.Module.Color==mod);
  cat(' ',length(rowind),'probes in module\n');
  submatrix = geneInfo[rowind,];
  orderind = order(submatrix[,colind],decreasing=TRUE);
  if (length(rowind) < numgenesingraph) {
    numgenesingraph = length(rowind);
    numconnections2keep = numgenesingraph * (numgenesingraph - 1);
  }
  cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
  submatrix = submatrix[orderind[1:numgenesingraph],];
  matchind = match(submatrix$GeneSymbol,colnames(datExpr));
  reducedTOM = TOM.matrix[matchind,matchind];
  
  orderind = order(reducedTOM,decreasing=TRUE);
  connections2keep = orderind[1:numconnections2keep];
  reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
  reducedTOM[connections2keep] = 1;
  
  gA <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
  gB <- graph.adjacency(as.matrix(reducedTOM[11:25,11:25]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutCircle <- rbind(layout.circle(gA)/2,layout.circle(gB))
  
  g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
  
  plot(g1,
       edge.color=adjustcolor(mod, alpha.f=0.25),
       edge.alpha=0.25,
       vertex.color=adjustcolor(mod, alpha.f=0.75),
       vertex.label=as.character(submatrix$GeneSymbol),
       vertex.label.cex=2.2,
       vertex.label.dist=1.1,
       vertex.label.degree=-pi/4,
       vertex.label.color="black",
       #vertex.frame.color='black',
       layout= jitter(layoutCircle),
       vertex.size=6,
       main=paste(mod,"module")
  )
}
# dev.off()

save(multiExpr, softPower, nSets, setLabels, shortLabels, file = "RData/Beta_multiExpre_Network.RData")
write.csv(reducedTOM,file=paste0('output/reduced_matrix.csv'))


# # next: further analysis
