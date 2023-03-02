## SET UP WORKING ENVIRONMENT ##
setwd('YOURWORKDIR')
rm(list = ls()); gc()
library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(BiocParallel)
library(infercnv)

## STEP 01: load in files and prepare data for inferCNV pipeline ##

pbmc = readRDS(file = 'pbmc.rds')
k20 = readRDS(file = 'K20.rds')
pdNormal = k20$pdNormal; pdmalignant = k20$pdMalignant; pdIntermediate = k20$pdIntermediate
counts = GetAssayData(pbmc, slot = 'counts')
meta = pbmc@meta.data
meta$barcode = rownames(meta)
metause = meta[, c('barcode', 'predicted')]
view(metause)
metause$predicted = 'filtered'
metause[which(metause$barcode %in% pdNormal),]$predicted = 'predicted non-neoplastic cells'
metause[which(metause$barcode %in% pdmalignant),]$predicted = 'predicted neoplastic cells'
metause[which(metause$barcode %in% pdIntermediate),]$predicted = 'predicted intermediate cells'

geneorder = readRDS(file = './gencode.vM31.pos.rds')
geneorder=geneorder[!duplicated(geneorder$gene_name),]
table(duplicated(geneorder$gene_name))
rownames(geneorder) = geneorder$gene_name 

genecommon = intersect(rownames(geneorder),rownames(counts))
counts = counts[genecommon,]
geneorder = geneorder[genecommon,]

geneorder=geneorder[, c('seqnames', 'start', 'end')]; names(geneorder) = c('chr', 'start', 'stop')

counts = as.matrix(counts)
metause.tmp = metause[, c('predicted')]
metause.tmp = as.data.frame(metause.tmp)
rownames(metause.tmp) = rownames(metause)

## STEP 02: InferCNV pipeline ##

infercnv_obj<-CreateInfercnvObject(raw_counts_matrix = counts,
                                   annotations_file = metause.tmp,
                                   gene_order_file = geneorder,
                                   ref_group_names = c('predicted non-neoplastic cells'))
options(scipen = 100)
options(expressions = 500000)
infercnv_obj<-run(infercnv_obj = infercnv_obj,
                  cutoff = 0.1, # cutoff=1 works well for smart-seq2, cutoff=0.1 works well for 10x
                  out_dir = './result_0215_refined',
                  cluster_by_groups = T,
                  denoise = T,
                  HMM = T,
                  analysis_mode = 'sample',
                  num_threads = 12,
                  num_ref_groups = NULL,
                  output_format = 'pdf',
                  scale_data = F,
                  write_expr_matrix = T)