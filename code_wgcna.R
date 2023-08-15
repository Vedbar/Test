# This program is for WGCNA on snoRNA
# Data from normalized Proteomic data and DESeq2 normalized snoRNA are merged

rm(list=ls())
library("WGCNA")
setwd("~/LTS/John_WGCNA/Jijiwa")
dir()
library("flashClust")
library("doParallel")
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

data= read.csv("snoRNA_protein_wcgna.csv", header = TRUE, row.names = 1)
head(data)
rowSums(data)
dim(data) #581 X 20
# Include rows only if non-empty cells is above 10 out of 20 cells
data_filtered <- data[rowSums(data != 0) >= 10, ] 
dim(data_filtered) # 475 X 20
head(data_filtered)
################################################################

gene.names=rownames(data_filtered)
Tdata = t(data_filtered)

#powers = c(c(1:10), seq(from = 12, to=20, by=2));
powers = c(1:11);
sft=pickSoftThreshold(Tdata,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 5 ;
#calclute the adjacency matrix
adj= adjacency(Tdata,type = "signed", power = softPower);
TOM=TOMsimilarityFromExpr(Tdata,networkType = "signed", TOMType = "signed", power = softPower);
colnames(TOM) =rownames(TOM) =gene.names
dissTOM=1-TOM

# Module detection
#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");
#plot the resulting clustering tree (dendrogram)
par(mfrow = c(1,1));
plot(geneTree, xlab="", sub="",cex=0.3);

# Set the minimum module size
minModuleSize = 15;

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = TRUE, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#set the diagonal of the dissimilarity to NA 
diag(TOM) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
sizeGrWindow(7,7)
TOMplot(TOM, geneTree, as.character(dynamicColors[dynamicColors]))

# Extract modules
module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=gene.names[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}


#####################################################################################
#####################################################################################
#####################################################################################
# MEList = moduleEigengenes(Tdata, colors = dynamicColors)
# MEs = MEList$eigengenes

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(Tdata, colors = dynamicColors)$eigengenes
# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
library(magrittr)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)
MEs0$treatment <- factor(MEs0$treatment, levels = c("HB1","HB2","HB3","HB4","HB6","HB7","HB10","HB11","HB12","HB13","HB14",
                           "CA1","CA2","CA3","CA4","CA5","CA6","CA7","CA8","CA9"))

library(dplyr)
library(tidyr)
# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

######################################################
# Alternative way
# mME = MEs0 %>%
#   pivot_longer(-treatment) %>%
#   mutate(
#     name = gsub("ME", "", name),
#     name = factor(name, levels = module_order)
#   )

MEList = moduleEigengenes(Tdata, colors = dynamicColors)
MEs = MEList$eigengenes
colnames(MEs)
#colnames(MEs)= gsub("ME", "", colnames(MEs))
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))


######################################################

t(MEs)
MEList$eigengenes
MEList$validColors

# Example data (replace this with your actual data)
#data_matrix <- t(MEs_modified)
data_matrix <- t(MEs)
#module_colors <- rep(c("red", "blue"), each = 5)

#sample_names <- paste("Sample", 1:10)
sample_names <- colnames(data_matrix)

# Alternative way
# mME = MEs0 %>%
#   pivot_longer(-treatment) %>%
#   mutate(
#     name = gsub("ME", "", name),
#     name = factor(name, levels = module_order)
#   )

library("gplots")
# Plot the heatmap
# Move the heatmap plot up by setting custom margins using the layout function
heatmap.2(data_matrix,
          Rowv = TRUE,   # Disable row clustering
          Colv = FALSE,   # Disable column clustering
          trace = "none", # Remove trace lines
          dendrogram = "row", # Remove dendrograms
          #labRow = module_colors, # Row labels are module colors
          labCol = sample_names,  # Column labels are sample names
          cexRow = 1, # Adjust the size of the row labels
          cexCol = 1.2, # Adjust the size of the column labels
          col = colorRampPalette(c("blue", "white", "red"))(100),  # Choose your desired color scale
          # main = "Module Colors and Samples",
          
          # Row/Column Labeling
          labRow = gsub("ME", "", rownames(data_matrix)),
          margins = c(7, 7),
          adjRow = c(0,NA),
          adjCol = c(NA,0),
          offsetRow  = 0.1,
          offsetCol = 1.2,
          srtCol=45,
          key = TRUE,
          keysize = 1.5,
          density.info = "none",
          key.title = "corr",
          key.ylab = "",
          key.xlab = "")

###################################################################################
#####################################################################################
#####################################################################################

library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(org.Hs.eg.db)
head( keys(org.Hs.eg.db, keytype="REFSEQ") ) # e.g, "NM_130786", "NP_570602"
head( keys(org.Hs.eg.db, keytype="GENENAME") ) # "alpha-1-B glycoprotein", "alpha-2-macroglobulin"
head( keys(org.Hs.eg.db, keytype="SYMBOL") ) # "A1BG", "A2M"
head( keys(org.Hs.eg.db, keytype="ENSEMBL") ) # "ENSG00000121410" "ENSG00000175899" 
columns(org.Hs.eg.db) ## display the columns
keytypes(org.Hs.eg.db) ## list supported key types

gene_list= gene.names[which(dynamicColors=="black")] # "MASP2", "SPTA1" 
# Perform GO enrichment analysis and save the results
go_enrichment <- enrichGO(gene      = gene_list,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",
                          ont          =  "ALL", # e.g. "BP", "MF", "CC"
                          pAdjustMethod = "BH",
                         # readable = TRUE,
                         # pvalueCutoff = 0.05,
                          qvalueCutoff = 1e-05)

go_enrichment
go_black <- data.frame(go_enrichment)
go_black = go_black[order(go_black$pvalue),]
head(go_black)
write.csv(go_black,"GO_black.csv")
dotplot(go_enrichment,
        #orderBy = "p.adjust", 
        #x= 'Count',
        showCategory=10)

emapplot(go_enrichment, showCategory = 2)
gene_list
ids <-bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db, drop = TRUE)
# remove duplicate
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
head(dedup_ids)
head(data)
dedup_ids$Y = dedup_ids$ENTREZID

data[rownames(data) %in% dedup_ids$SYMBOL, ][1:11]
data[rownames(data) %in% dedup_ids$SYMBOL, ][12:20]

gene_sum <- aggregate(. ~ Gene, data = selected_data, sum)



# MEs0 <- moduleEigengenes(Tdata, colors = dynamicColors)$eigengenes
# Tdata
kk2 <- gseKEGG(geneList     = dedup_ids,
               organism     = "hsa",
               nPerm        = 100,
               minGSSize    = 3,
               maxGSSize    = 10,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

##########################################################################
for (color in module_colors){
module=gene.names[which(dynamicColors==color)]
gene_list= gene.names[which(dynamicColors==color)] # "MASP2", "SPTA1" 
# Perform GO enrichment analysis and save the results
go_enrichment <- enrichGO(gene      = gene_list,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",
                          ont          =  "ALL", # e.g. "BP", "MF", "CC"
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1e-05)
go_data <- data.frame(go_enrichment)
write.csv(go_data,paste("GO_",color, ".csv",sep=""))
}

##########################################################################

OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene
## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(go_enrichment, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange = NULL, 
         vertex.label.font=6)

cnetplot(go_enrichment)

#############################
data(geneList, package='DOSE')
de2 <- names(geneList)[1:2]
yy2 <- enrichKEGG(de2)
cnetplot(yy2)


module_enrichment_results <- lapply(gene.names, function(gene_list) {
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(gene      = gene_list,
                            OrgDb        = org.Hs.eg.db,
                            keyType      = "SYMBOL",
                            ont          = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.02)
  
  # Perform KEGG pathway enrichment analysis
  kegg_enrichment <- enrichKEGG(gene         = gene_list,
                                organism    = 'hsa', # #aga' Organism code for Anopheles gambiae
                                keyType = "SYMBOL",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2)
  
  # Return both GO and KEGG enrichment results for the module
  return(list(GO = go_enrichment, KEGG = kegg_enrichment))
})

# Perform GO enrichment analysis using DOSE package
go_enrichment <- enrichGO(gene     = gene_list,
                          OrgDb    = org.Hs.eg.db,
                          # keyType  = "ENTREZID",
                          keyType = "SYMBOL", # "GENENAME",
                          ont      = "BP",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH")
go_enrichment
cluster_summary <- data.frame(go_enrichment)
return(go_enrichment)
})


# #####################################
# library(KEGGREST)
# library(KEGGprofile)
# gene_symbols <- c("TP53", "BRCA1", "EGFR")
# gene_list
# # Convert gene symbols to KEGG gene IDs
# converted_ids <- keggConv(gene_symbols)
# # Print the result
# print(converted_ids)
# #####################################