
#1. Enrichment analysis
# In this pipeline we use clusterProfiler

## package
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(carData)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)


# load data --------------------------------------------------------------------
# setwd("F:/$%/@/GO")

load("Gamma_module_genes_2.RData")
# write.csv(gene_Mod_blue, file = "results\\gene_M1Blue.csv", row.names = F)
# write.csv(gene_Mod_brown, file = "results\\gene_M1Brown.csv", row.names = F)
# write.csv(gene_Mod_grey, file = "results\\gene_M1Gray.csv", row.names = F)
# write.csv(gene_mod_turquoise, file = "results\\gene_M1Gray.csv", row.names = F)

# color data --------------------------------------------------------------

# gene = Alpha_Mod_blue
# gene = gene_Mod_blown
gene = Gamma_Mod_grey
# gene <- delta_Mod_turq



# genes symbol ----------------------------------------------------------
genelist <- as.character(gene)
eg <- bitr(genelist, 
           fromType="SYMBOL", 
           toType=c("ENTREZID","GENENAME"), 
           OrgDb="org.Hs.eg.db"); 
head(eg)

geneList <- eg$ENTREZID

# go ----------------------------------------------------------------------

go <- enrichGO(gene = geneList,
               OrgDb = org.Hs.eg.db,
               ont='ALL',            # replace with BP, CC and MF
               pAdjustMethod = 'BH',
               ## Blue
               # pvalueCutoff = 0.2,
               # qvalueCutoff = 0.5,
               ## Brown
               # pvalueCutoff = 0.05,
               # qvalueCutoff = 0.05,
               # Grey
               pvalueCutoff = 0.2,
               qvalueCutoff = 0.5,
               # # Tur
               # pvalueCutoff = 0.05,
               # qvalueCutoff = 0.05,
               keyType = 'ENTREZID')
head(go)
# write.csv(go@result, file = "res_delta\\go_MBlue.csv", row.names = F) 
# write.csv(go@result, file = "res_delta\\go_Mbrown.csv", row.names = F) 
write.csv(go@result, file = "res_gamma\\go_MGrey.csv", row.names = F)
# write.csv(go@result, file = "res_delta\\go_MTur.csv", row.names = F)

# View(go@result)
# visualization ----------------------------------------------------------------
barplot(go, showCategory=10,drop=T, title=NULL)    # showCategory=10,title="Enrichment GO"
dotplot(go, showCategory=10) 



## 2022.6.5 add
enrichKEGG(genelist, 
           organism = 'hsa',
           keyType = 'kegg', 
           # organism = 'human',
           pvalueCutoff = 0.05, 
           pAdjustMethod = 'BH', 
           minGSSize = 3, 
           maxGSSize = 500, 
           qvalueCutoff = 0.2, 
           use_internal_data = FALSE)


barplot(enrichKK,showCategory=20)
dotplot(enrichKK)

#pathway 
browseKEGG(kegg, "hsa00770") 

# 2. Differential Co-expression network analysis
#plot differential coexpression Gene network
# We used MODA package to perform thi step

# BiocManager::install("MODA")
library(MODA)

load("Beta_datExpre.RData")   #replace by other cell types data

ResultFolder = 'Beta' # where middle files are stored
CuttingCriterion = 'Density' # could be Density or Modularity
indicator1 = 'X'     # indicator for data profile 1
indicator2 = 'Y'      # indicator for data profile 2
specificTheta = 0.1 #threshold to define condition specific modules
conservedTheta = 0.1#threshold to define conserved modules
##modules detection for network 1
BetaModules <- WeightedModulePartitionHierarchical(datExpr,
                                                   ResultFolder,
                                                   indicator1,
                                                   CuttingCriterion)
# png('Beta/heatdatExpr.png')
tiff("Beta/Beta hist genes per cell.tiff", units="in", width=6, height=4, res=300)
heatmap(cor(as.matrix(datExpr)))
dev.off()
