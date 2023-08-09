rm(list = ls())
gc()
setwd("YOURWORKDIR")
library(tidyverse)
library(data.table)
library(stringr)
library(dplyr)
library(magrittr)
source("./copyscat_refined_functions.R")
library(CopyscAT)
library(BSgenome.Mmusculus.UCSC.mm10)
library(scales)

## initialize the environment ##
initialiseEnvironment(genomeFile = "./mm10_chrom_sizes.tsv",
                      cytobandFile = "./mm10_1e+6_cytoband_densities_granges.tsv",
                      cpgFile = "./mm10_1e+6_cpg_densities.tsv",
                      binSize = 1e+6,
                      minFrags = 1e4,
                      cellSuffix = c("-1"),
                      lowerTrim = 0.3,
                      upperTrim = 0.9)

#SET OUTPUT DEFAULT DIRECTORY AND NAME
setOutputFile("./100wbin/", "lung_cnv_total_1e+06")

#PART 1: INITIAL DATA NORMALIZATION
#step 1 normalize the matrix
counts <- readRDS(file = "./total_pseudocell_100wbin.rds")
counts <- as.data.frame(t(counts))
scData <- counts
rownames(scData) <- paste0(rownames(scData), '-1')

scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=500,upperFilterQuantile = 1, dividingFactor = 1)
summaryFunction <- cutAverage

scData_collapse <- collapseChrom3N(scData_k_norm, summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 100,logTrans = FALSE,tssEnrich = 4,logBase=2, minCPG=300,powVal=0.73)
#apply additional filters
scData_collapse_temp <- filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.1)
scData_collapse <- scData_collapse_temp

######## NMF #######
if(F)
{
  inputMatrix <- scData_collapse
  nmfComponents <- 12
  methodHclust <- 'ward.D2'
  cutClusters <- 4
  outputHeatmap <- TRUE
  #uses NMF and fast hclust packages
  message("Running NMF")
  res <- nmf(column_to_rownames(inputMatrix,var="chrom"), c(nmfComponents), "brunet", seed="nndsvd",.stop=nmf.stop.threshold(0.1),maxIter=2500)
  message("Computing clusters")
  tscale <- scale(x = t(coef(res)), center = TRUE)
  dist <- dist(tscale, method = "euclidean")
  clusts <- fastcluster::hclust(dist, method = methodHclust)
  cell_assigns <- cutree(clusts, k = cutClusters)
  table(cell_assigns)
  plot(clusts)
  if (outputHeatmap == TRUE)
  {
    pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_heatmap.pdf"),width=6,height=6)
    heatmap.2(coef(res),Rowv=FALSE,Colv=as.dendrogram(clusts),dendrogram="column",density.info="none",trace="none",scale="none",labCol=FALSE,col=colorRampPalette(viridis(5)),symkey=FALSE,useRaster=TRUE,ColSideColors=viridis(n=length(unique(as.character(cell_assigns))))[cell_assigns])
    legend("topright",fill=viridis(n=3),x.intersp = 0.8,y.intersp=0.8,legend=unique(as.character(cell_assigns)),horiz = FALSE,cex = 0.9,border=TRUE,bty="n")
    dev.off()
  }
  by_cluster <- left_join(rownames_to_column(data.frame(t(column_to_rownames(inputMatrix,"chrom"))),var="barcode"),
                        rownames_to_column(data.frame(cluster=cell_assigns),var="barcode"),by="barcode")
  #now estimate the residuals etc for each cluster
  if (outputHeatmap == TRUE)
  {
    pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_violinplot.pdf"),width=12,height=6)
    print(ggplot(by_cluster %>% gather(starts_with("chr"),key="chrom",value="value"),aes(chrom,value)) + geom_violin() + theme(axis.text.x=element_text(angle=-90,vjust = 0.5,hjust = 0, color = "#000000",size = 6)) + facet_wrap(~cluster))
    dev.off()
  }
  ### assignment ###
  cell.df <- as.data.frame(cell_assigns)
  cell.df$group <- ""
  cell.df$group[grep("adj", rownames(cell.df))] <- "adj"
  cell.df$group[grep("tumor", rownames(cell.df))] <- "tumor"
  cell.df$group[grep("wt", rownames(cell.df))] <- "wt"
  df <- as.data.frame(table(cell.df$cell_assigns, cell.df$group))
  wt.prop <- df[df$Var2 == "wt", ]$Freq / sum(df[df$Var2 == "wt",]$Freq)
  names(wt.prop) <- unique(df$Var1)
  target_cluster_1 <- names(wt.prop)[which(wt.prop >= 0.8)]
  print("potential normal hmf cluster:")
  print(target_cluster_1)
  target_cluster_1 <- data.frame(cluster=as.numeric(target_cluster_1))
  cell.df[which(cell.df$cell_assigns != target_cluster_1$cluster & cell.df$group == 'wt'),]$cell_assigns = target_cluster_1$cluster
  cell_assigns <- cell.df$cell_assigns
  names(cell_assigns) <- rownames(cell.df)
  tumor_cell_ids <- names(which(cell_assigns != target_cluster_1))
  normalCluster <- unlist(target_cluster_1)
  normal_cell_ids <- names(which(cell_assigns == unlist(normalCluster)))
  nmf_results <- list(cellAssigns = cell_assigns,normalBarcodes = normal_cell_ids, clusterNormal = normalCluster)
  write.table(x = rownames_to_column(data.frame(nmf_results$cellAssigns), var="Barcode"),file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_clusters.csv"), quote=FALSE, row.names = FALSE, sep=",")
  print(paste("Normal cluster is: ",nmf_results$clusterNormal))
}

median_iqr <- computeCenters(scData_collapse %>% dplyr::select(chrom,nmf_results$normalBarcodes),summaryFunction=summaryFunction)
candidate_cnvs <- identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = T,minMix=0.001,deltaMean = 0.1,deltaBIC2 = 0.5,bicMinimum = 0.01, uncertaintyCutoff = 0.55,summaryFunction=summaryFunction,maxClust = 6,mergeCutoff = 3,IQRCutoff =0.1,medianQuantileCutoff = -1,normalCells=nmf_results$normalBarcodes, subsetSize = 40) 
candidate_cnvs_clean <- clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.5, maxClust = 6)

final_cnv_list <- annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "_clean_cnv",sdCNV = 0.5,filterResults=FALSE, minAlteredCells = 5)
print("potentially varied CNV:")
print(dim(final_cnv_list[[3]])[2])
