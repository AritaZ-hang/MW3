## SET UP WORKING ENVIRONMENT ##

rm(list = ls())
library(ArchR)
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(data.table)
options(scipen = 200)
library(FNN)

setwd('YOURWORKDIR')

## STEP 01: Build imputed matrix for scATAC-seq data ##
## We've done this by intergating it with scRNA-seq data using Seurat's FindTransferAnchors function
oligo_atac = readRDS('oligo_atac.rds')
oligo_rna = readRDS('oligo_rna.rds')
oligo_cds = readRDS('oligo_cds.rds') # obtained by Monocle2
rna_counts = GetAssayData(oligo_rna, slot = 'counts')

oligo_geneScore = getMatrixFromProject(oligo_atac, useMatrix = 'GeneScoreMatrix')
oligo_geneMat = oligo_geneScore@assays@data@listData[["GeneScoreMatrix"]]
rownames(oligo_geneMat) = oligo_geneScore@elementMetadata@listData[["name"]]

genes_intersect = intersect(rownames(oligo_geneMat), rownames(oligo_rna))

oligo_geneMat = oligo_geneMat[genes_intersect,]
oligo_rna = oligo_rna[genes_intersect,]

######################### Deal with gene scores ###################################
if(T)
{
  oligo_seurat = CreateSeuratObject(counts = oligo_geneMat, min.cells = 3,project = "oligo")
  oligo_seurat@assays$RNA@data = oligo_seurat@assays$RNA@counts
  oligo_seurat = NormalizeData(object = oligo_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  oligo_seurat = FindVariableFeatures(object = oligo_seurat, mean.function = ExpMean, dispersion.function = LogVMR)
  hv.genes.atac= head(VariableFeatures(oligo_seurat),2000)
  oligo_seurat = ScaleData(object = oligo_seurat,features = hv.genes.atac, vars.to.regress = c("nCount_RNA"))
  oligo_seurat = RunPCA(object = oligo_seurat, pc.genes = hv.genes.atac, pcs.compute = 100, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
  ElbowPlot(object =oligo_seurat, ndims = 50)
  
  oligo_seurat = FindNeighbors(oligo_seurat,dims=1:50)
  oligo_seurat = FindClusters(oligo_seurat, resolution = 1)
  oligo_seurat = RunUMAP(object = oligo_seurat, reduction = "pca", dims= 1:50, reduction.name = "umap") 
}

all.hv_genes = unique(c(VariableFeatures(oligo_seurat), VariableFeatures(oligo_rna)))

transfer.anchors.oligo = FindTransferAnchors(reference = oligo_rna, query = oligo_seurat, reference.assay = 'RNA', query.assay = 'RNA',reduction = 'cca', k.anchor = 30, features = all.hv_genes, npcs = 50)

## STEP 02: Mapping ATAC cells with neighbor RNA cells in the common CCA space by FNN ## 

ImputeCoembed = function(rna_object, atac_object, transfer.anchors)
{
  genes.use = transfer.anchors@anchor.features
  refdata = GetAssayData(rna_object, assay = 'RNA', slot = 'data')[genes.use,]
  imputation = TransferData(transfer.anchors, refdata, weight.reduction = 'cca', l2.norm = T)
  refdata.matrix = t(as.matrix(refdata))
  imputation.matrix = GetAssayData(imputation, slot = 'data')
  imputation.matrix = t(as.matrix(imputation.matrix))
  
  return(list(refdata.matrix, imputation.matrix))
  
}
sampleNearestNeighbors = function(cell, knn_object, N = 30, output = 'indices')
{
  if(output == 'indices')
  {
    cells = as.integer(as.data.frame(knn_object[['nn.index']])[cell,])
  }
  else
  {
    cells = as.integer(as.data.frame(knn_object[['nn.cells']])[cell,])
  }
  dists = as.double(as.data.frame(knn_object[['nn.dist']])[cell,])
  prob_fun = eval(parse(text = "(1/exp(dists)) * (1/sum(1/exp(dists)))"))
  sample(x = cells,
         size = N,
         replace = F,
         prob = prob_fun)
}
FnnMatching = function(impute.ref, kn = 50, random = 20)
{
  knn.grouping = get.knnx(data = impute.ref[[1]],
                          query = impute.ref[[2]],
                          algorithm = 'kd_tree', 
                          k = kn)
  temp.index = c()
  for(i in 1:length(rownames(impute.ref[[2]])))
  {
    idx = sampleNearestNeighbors(cell = i, knn_object = knn.grouping, N = random)
    temp.index = rbind(temp.index, idx)
  }
  temp.index = as.data.frame(temp.index)
  rownames(temp.index) = rownames(impute.ref[[2]])
  rna_index.df = data.frame(rna_cells = rownames(impute.ref[[1]]),
                            index = 1:length(rownames(impute.ref[[1]])))
  final.index = c()
  for(i in 1:random)
  {
    atac_index.df = data.frame(atac_cells = rownames(temp.index),
                               index = temp.index[,i])
    t = dplyr::inner_join(atac_index.df, rna_index.df, by = 'index')[, 3]
    final.index = cbind(final.index, t)
  }
  final.index = as.data.frame(final.index)
  rownames(final.index) = rownames(impute.ref[[2]])
  names(final.index) = 1:random
  return(final.index)
}
oligo.impute.ref = ImputeCoembed(oligo_rna, oligo_seurat, transfer.anchors.oligo)
oligo_matching_df = FnnMatching(oligo.impute.ref, 100, 50)
oligo_matching_df$atac_cells = rownames(oligo_matching_df)

## STEP 03: Compute average Pseudotime ##
oligo_rna_pseudotime = oligo_cds@phenoData@data

MeanPseudotime = function(matching_df, rna_pseudotime, random = 20)
{
  atac_cells = unique(matching_df$atac_cells)
  atac_pseudotime = c()
  for(i in atac_cells)
  {
    subset = matching_df[matching_df$atac_cells == i,]
    mean_pseudo = mean(rna_pseudotime[match(subset[i,1:random], rownames(rna_pseudotime)),]$Pseudotime_regress)
    atac_pseudotime = c(atac_pseudotime, mean_pseudo)
  }
  atac_pseudotime = as.data.frame(atac_pseudotime)
  rownames(atac_pseudotime) = atac_cells
  return(atac_pseudotime)
}

oligo_atac_pseudotime = MeanPseudotime(oligo_matching_df, oligo_rna_pseudotime, random = 50)

## STEP 03: Build corresponding pseudobulks for scRNA-seq data ## 
BuildPseudobulk = function(matching_df, counts)
{
  sum = c()
  for(i in 1:dim(matching_df)[1])
  {
    cells_to_use = as.character(matching_df[i, 1:50])
    tmp_counts = counts[, cells_to_use]
    tmp = as.vector(Matrix::rowSums(tmp_counts))
    sum = cbind(sum, tmp)
  }
  sum = as.data.frame(sum)
  rownames(sum) = rownames(counts)
  names(sum) = paste0('Pseudobulk_', 1:dim(matching_df)[1])
  return(sum)
}

rna_pseudobulk = BuildPseudobulk(new_match_df, rna_counts)
rna_pseudobulk_pseudotime = MeanPseudotime(new_match_df, oligo_rna_pseudotime, random = 50)
names(rna_pseudobulk_pseudotime) = 'pseudotime'
rownames(rna_pseudobulk_pseudotime) = paste0('Pseudobulk_', 1:length(cells_needed))

save(rna_pseudobulk, rna_pseudobulk_pseudotime, file = 'oligo_rna_pseudobulk.RData')
saveRDS(oligo_atac_pseudotime, file = 'oligo_fnn_pseudotime.rds')
