## SET UP WORKING ENVIRONMENT ##
rm(list=ls()); gc()
setwd('YOURWORKDIR')
library(tidyverse)
library(data.table)
library(stringr)
library(ArchR)
addArchRGenome('mm10')
library(dplyr)
library(magrittr)
source('CopyscAT_refinedFunctions.R')

library(CopyscAT)
library(BSgenome.Mmusculus.UCSC.mm10)

## STEP 01: initialize the environment ##

## initialize the environment ##
initialiseEnvironment(genomeFile = './mm10_chrom_sizes.tsv',
                      cytobandFile = './mm10_5e+06_cytoband_densities_granges.tsv',
                      cpgFile = './mm10_5e+06_cpg_densities.tsv',
                      binSize = 5e+06,
                      minFrags = 1e4,
                      cellSuffix = c('-1'),
                      lowerTrim=0.4,
                      upperTrim = 0.9)
setOutputFile("./","lung_cnv_total_5e+06")

## STEP 02: Normalize the matrix ##

scData = readInputTable("./lung.wt.tumor.sort.500w.tsv")
scData_k_norm =  normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=125,upperFilterQuantile = 0.95, dividingFactor = 1)
summaryFunction = cutAverage
scData_collapse = collapseChrom3N(scData_k_norm, summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 200,logTrans = FALSE,tssEnrich = 4,logBase=2, minCPG=300,powVal=0.73)
scData_collapse_temp = filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.1)
scData_collapse = scData_collapse_temp

## STEP 03: using NMF to distinguish non-neoplastic cells
if(T)
{
  #uses NMF and fast hclust packages
  message("Running NMF")
  res = nmf(column_to_rownames(scData_collapse,var="chrom"), 12, "lee", seed="nndsvd",.stop=nmf.stop.threshold(0.1),maxIter=2500)
  message("Computing clusters")
  tscale = scale(x=t(coef(res)),center=TRUE)
  dist = dist(tscale,method="euclidean")
  
  clusts = fastcluster::hclust(dist,method = 'ward.D')
  cell_assigns = cutree(clusts,h=max(clusts$height)*cutHeight)
  table(cell_assigns)
  
  cell.df = as.data.frame(cell_assigns)
  cell.df$group = str_sub(rownames(cell.df), start = 1, end = 3)
  
  ## Add Tissue info
  RT = fread(file = 'RT.csv', header = F)$V1
  RT.ls = data.frame(RT=RT, tissue='')
  RT.ls[289:336,]$tissue='Tumor'
  RT.ls[337:384,]$tissue='Adj'
  
  barcodes.ls = data.frame(cells = rownames(cell.df[which(cell.df$group == 'COL'),]),
                           RT = str_sub(rownames(cell.df[which(cell.df$group == 'COL'),]), start = -12, end = -3))
  
  barcodes.ls.ref=dplyr::inner_join(barcodes.ls,RT.ls,by='RT')
  table(barcodes.ls.ref$tissue)
  
  cell.df[barcodes.ls.ref$cells[which(barcodes.ls.ref$tissue == 'Tumor')],]$group = 'COL_tumor'
  cell.df[barcodes.ls.ref$cells[which(barcodes.ls.ref$tissue == 'Adj')],]$group = 'COL_adj'
  
  df = as.data.frame(table(cell.df$cell_assigns, cell.df$group))
  bat.prop = df[df$Var2 == 'bat',]$Freq / sum(df[df$Var2 == 'bat',]$Freq)
  names(bat.prop) = unique(df$Var1)
  bat.prop
  
  target_cluster_1 = names(bat.prop)[which(bat.prop >= 0.8)]
  print(target_cluster_1)

  target_cluster_1 = data.frame(cluster=as.numeric(target_cluster_1))
  cell.df[which(cell.df$cell_assigns != target_cluster_1$cluster & cell.df$group == 'bat'),]$cell_assigns = target_cluster_1$cluster
  cell_assigns = cell.df$cell_assigns
  names(cell_assigns) = rownames(cell.df)
  table(cell_assigns)
  
  tumor_cell_ids=names(which(cell_assigns!=target_cluster_1))
  normalCluster=unlist(target_cluster_1)
  normal_cell_ids=names(which(cell_assigns==unlist(normalCluster)))
  
  nmf_results = list(cellAssigns=cell_assigns,normalBarcodes=normal_cell_ids,clusterNormal=normalCluster)
  
}

write.table(x=rownames_to_column(data.frame(nmf_results$cellAssigns),var="Barcode"),file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_clusters.csv"),quote=FALSE,row.names = FALSE,sep=",")
print(paste("Normal cluster is: ",nmf_results$clusterNormal))
table(nmf_results$cellAssigns)

## STEP 04: Compute copy number alterations
  
median_iqr = computeCenters(scData_collapse %>% dplyr::select(chrom,nmf_results$normalBarcodes),summaryFunction=summaryFunction)
candidate_cnvs = identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,propDummy=0.4,minMix=0.01,deltaMean = 0.03,deltaBIC2 = 0.25,bicMinimum = 0.1, subsetSize=800,fakeCellSD = 0.15, uncertaintyCutoff = 0.65,summaryFunction=summaryFunction,maxClust = 4,mergeCutoff = 3,IQRCutoff = 0.25,medianQuantileCutoff = -1,normalCells=nmf_results$normalBarcodes) 
candidate_cnvs_clean = clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff= 0.8) #= 1.5)
final_cnv_list = annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "_clean_cnv",sdCNV = 0.9,filterResults=FALSE)

