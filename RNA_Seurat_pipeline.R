## SET UP WORKING ENVIRONMENT ##
setwd('YOURWORKDIR')
set.seed(1)
library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(scDblFinder)
library(tidyverse)
library(DESeq2)
library(BiocParallel)

rm(list=ls()); gc()

## STEP 01: read dge file from the folder ##
dge = fread("lung_rmbatch.csv")
dge = as.data.frame(dge)
rownames(dge) = dge$V1
dge = dge[,-1]
colnames(dge) = paste0("TISSUENAME",".",colnames(dge))

## STEP 02: Create Seurat Object ##
pbmc = CreateSeuratObject(counts = dge, min.cells = 3, project = "lung")

## STEP 03: Quality Control ##
summary(colSums(dge))
hist(colSums(dge),breaks = 100)
abline(v=500)
hist(colSums(dge>0),breaks = 100)
abline(v=200)
summary(colSums(dge>0))
mito.genes = grep(pattern = "^mt-", x = rownames(x = pbmc@assays$RNA), value = TRUE)
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^mt-")

pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10
               &  nCount_RNA < 20000)

## STEP 04: Preprocessing ##
pbmc = FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR)
hv.genes = head(VariableFeatures(pbmc),2000)
a = setdiff(VariableFeatures(pbmc),mito.genes) ## remove mt.genes
hv.genes = a

pbmc = ScaleData(object =pbmc,features = hv.genes, vars.to.regress = c("percent.mt","nCount_RNA"))
pbmc = RunPCA(object =pbmc, pc.genes = hv.genes, pcs.compute = 100, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
ElbowPlot(object =pbmc, ndims = 50)

## STEP 05: Find clusters and run UMAP ## 
pbmc = FindNeighbors(pbmc,dims=1:50)
pbmc = FindClusters(pbmc, resolution = 1)
pbmc = RunUMAP(object = pbmc, reduction = "pca", dims= 1:50, reduction.name = "umap") 

## STEP 06: Obtain markers for each clusters ## 
pbmc.markers = FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.1, logfc.threshold  = 0.25)