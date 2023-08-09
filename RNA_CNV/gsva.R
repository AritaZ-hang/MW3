# GSVA in pseudobulk versions
## set up working environment..
rm(list = ls())
gc()
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
library(org.Mm.eg.db)

setwd("YOURWORKDIR")
load("./adj_pseudobulks.RData")
load("./ca_pseudobulks.RData")
k50 <- readRDS(file = "./hclust_pseudobulks_k50_iterativeHclust_define.rds")

k50_predict <- data.frame(pseudobulks_name = c(k50$pdNormal, k50$pdMalignant, k50$pdIntermediate), 
                         predict_new = c(rep("pdNormal", length(k50$pdNormal)), rep("pdMalignant", length(k50$pdMalignant)), rep("pdIntermediate", length(k50$pdIntermediate))))
rownames(k50_predict) <- k50_predict$pseudobulks_name

pseudobulk_merge <- readRDS(file = "./total_pseudocell_kn200.rds")
pseudobulk_merge_seob <- CreateSeuratObject(counts = pseudobulk_merge,
meta = k50_predict)

geneset <- msigdbr(species = "Mus musculus", category = "H")
geneset <- subset(geneset,
select = c("gs_name", "gene_symbol")) %>% as.data.frame()
geneset <- split(geneset$gene_symbol, geneset$gs_name)

col_seq <- c("pdNormal", "pdIntermediate", "pdMalignant")

expr <- AverageExpression(pseudobulk_merge_seob,
assays = "RNA", slot = "counts", group.by = "predict_new")[[1]]
expr <- expr[rowSums(expr) > 0, ]
gsva.res <- gsva(expr,
  geneset, method = "gsva", kcdf = "Poisson", parallel.sz = 5)
gsva.res <- gsva.res[, col_seq]
p_all <- pheatmap::pheatmap(gsva.res, show_colnames = T, scale = 'none', angle_col = '45',
                           color = colorRampPalette(c('#6898C5', 'white', '#C9676F'))(100),
                           border_color = NA, width = 10, height = 8,
                           cluster_cols = F
)
gsva.res <- gsva.res[p_all$tree_row$order,]
pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "none",
                   color = colorRampPalette(c('#6898C5', 'white', '#C9676F'))(100), 
                   border_color = "white", width = 10, height = 8,
                   filename = "./all_hallmark_in_pseudobulks.pdf", cluster_cols = F, cluster_rows = F, 
                   angle_col = 0)

