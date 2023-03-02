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
## We've done this by intergating it with scRNA-seq data using Signac
oligo_atac = readRDS('oligo_atac.rds')
oligo_rna = readRDS('oligo_rna.rds')
oligo_cds = readRDS('oligo_cds.rds') # obtained by Monocle2

oligo_geneScore = getMatrixFromProject(oligo_atac, useMatrix = 'GeneScoreMatrix')
oligo_geneMat = oligo_geneScore@assays@data@listData[["GeneScoreMatrix"]]
rownames(oligo_geneMat) = oligo_geneScore@elementMetadata@listData[["name"]]

## Standard workflow as Signac recommends.
assay = CreateChromatinAssay(counts = oligo_mat, sep = c('_', '_'))
oligo_seurat = CreateSeuratObject(counts = assay, assay = 'peaks',project = 'ATAC')

oligo_seurat = RunTFIDF(oligo_seurat)
oligo_seurat = FindTopFeatures(oligo_seurat, min.cutoff = 'q0')
oligo_seurat = RunSVD(object = oligo_seurat)
oligo_seurat = RunUMAP(object = oligo_seurat,reduction = 'lsi',dims = 2:30)
oligo_seurat = FindNeighbors(object = oligo_seurat,reduction = 'lsi',dims = 2:30)
oligo_seurat = FindClusters(object = oligo_seurat,algorithm = 3,resolution = 1,verbose = FALSE)

oligo_seurat[['activity']] = CreateAssayObject(counts = oligo_geneMat)
DefaultAssay(oligo_seurat) = 'activity'
oligo_seurat = ScaleData(oligo_seurat, features = rownames(oligo_seurat))

transfer.anchors.oligo = FindTransferAnchors(reference = oligo_rna, query = oligo_seurat, reference.assay = 'RNA', query.assay = 'activity',reduction = 'cca', k.anchor = 20, features = VariableFeatures(oligo_rna))

## STEP 02: Mapping ATAC cells with neighbor RNA cells in the common CCA space by FNN ## 

ImputeCoembed = function(rna_object, atac_object, transfer.anchors)
{
  genes.use = VariableFeatures(rna_object)
  refdata = GetAssayData(rna_object, assay = 'RNA', slot = 'data')[genes.use,]
  imputation = TransferData(transfer.anchors, refdata, weight.reduction = atac_object[['lsi']], dims = 2:30)
  atac_object[['RNA']] = imputation
  DefaultAssay(atac_object) = 'RNA'
  refdata.matrix = t(as.matrix(refdata))
  imputation.matrix = GetAssayData(imputation, slot = 'data')
  imputation.matrix = t(as.matrix(imputation.matrix))
  return(list(refdata.matrix, imputation.matrix))
  
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
oligo_matching_df = FnnMatching(oligo.impute.ref, 100, 100)
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

oligo_atac_pseudotime = MeanPseudotime(oligo_matching_df, oligo_rna_pseudotime, random = 100)

## STEP 03: sampling cells ##

fnx = ecdf(oligo_atac_pseudotime$atac_pseudotime)
t1 = fnx(oligo_atac_pseudotime$atac_pseudotime)
sub_cells = sample(rownames(oligo_atac_pseudotime), size = 100, replace = F, prob = 1/t1 * (1/sum(1/exp(t1))))
sample.atac = as.data.frame(oligo_atac_pseudotime[sub.cells,])
rownames(sample.atac) = test; names(sample.atac) = 'atac_pseudotime'
new_match_df = oligo_matching_df[sub_cells,]
rna_counts = GetAssayData(oligo_b, slot = 'counts')

## STEP 04: Build corresponding pseudobulks for scRNA-seq data ## 
BuildPseudobulk = function(matching_df, counts)
{
  sum = c()
  for(i in 1:dim(matching_df)[1])
  {
    cells_to_use = as.character(matching_df[i, 1:20])
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
rna_pseudobulk_pseudotime = MeanPseudotime(new_match_df, oligo_rna_pseudotime, random = 100)
names(rna_pseudobulk_pseudotime) = 'rna_pseudobulk_pseudotime'

save(sample.atac, rna_pseudobulk, rna_pseudobulk_pseudotime, file = 'oligo_rna_pseudobulk_100.RData')
saveRDS(oligo_atac_pseudotime, file = 'oligo_fnn_pseudotime_100.rds')

save.image(file = 'oligo_fnn_mapping_100.RData')
