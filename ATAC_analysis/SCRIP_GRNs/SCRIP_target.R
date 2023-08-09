## SET UP WORKING ENVIRONMENT ##
rm(list = ls()); gc()
library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(patchwork)
library(data.tree)
library(gridExtra)
library(rlist)
library(scales)
library(dendextend)
library(tidytree)
library(ape)
library(clusterProfiler)
library(org.Mm.eg.db)
library(igraph)
library(ggraph)
library(tidygraph)
library(data.table)
library(SeuratDisk)

setwd('YOURWORKDIR')

## STEP 01: load in file & remove low quality targets..

Nkx2_1.TR = Read10X('./Nkx2-1/', gene.column = 1)
delete = data.frame(apply(Nkx2_1.TR, 1, sum))
colnames(delete) = 'sum'
delete$name = rownames(delete)
delete = delete[order(delete[,1], decreasing = F),]
summary(delete$sum)
delete.tmp = delete[delete$sum < 300,] 
delete_target = rownames(delete.tmp)
Nkx2_1.TR.filtered = Nkx2_1.TR[setdiff(rownames(Nkx2_1.TR), delete_target),]

## STEP 02: Create Seurat Object & pre-processing

Nkx2_1.TR_seu = CreateSeuratObject(counts = (Nkx2_1.TR.filtered))
Nkx2_1.TR_seu = NormalizeData(Nkx2_1.TR_seu)
Nkx2_1.allGenes = rownames(Nkx2_1.TR_seu)
Nkx2_1.TR_seu = ScaleData(Nkx2_1.TR_seu, features = Nkx2_1.allGenes)

## STEP 03: Change the @active.ident to CopyscAT's predictions.
## This is skipped here.

## STEP 04: Find potential in neoplastic cells & non-neoplastic cells respectively.
## We'll do this by differential expression.

cluster_markers_non = FindMarkers(Nkx2_1.TR_seu, 
                                  ident.1 = 'non-neoplastic cells',
                                  ident.2 = 'neoplastic cells',
                                  logfc.threshold = 0.1,
                                  only.pos = T)
clusters_markers_neo = FindMarkers(Nkx2_1.TR_seu,
                                   ident.1 = 'neoplastic cells',
                                   ident.2 = 'non-neoplastic cells',
                                   logfc.threshold = 0.1,
                                   only.pos = T)

non_target.nkx2_1 = rownames(cluster_markers_non)[which(cluster_markers_non$avg_log2FC > 0.25 & cluster_markers_non$p_val_adj < 0.05)]
neo_target.nkx2_1 = rownames(clusters_markers_neo)[which(clusters_markers_neo$avg_log2FC > 0.25 & clusters_markers_neo$p_val_adj < 0.05)]

## remove histone genes
non_target.nkx2_1 = non_target.nkx2_1[-grep(pattern='H[1-4][a-z]+', non_target.nkx2_1)]
#neo_target.nkx2_1 = neo_target.nkx2_1[-grep(pattern='H[1-4][a-z]+', neo_target.nkx2_1)]

## STEP 05: GO:BP Enrichment analysis ##

eg_non = bitr(non_target.nkx2_1, fromType = 'SYMBOL', toType = c('ENTREZID', 'ENSEMBL', 'SYMBOL'), OrgDb = org.Mm.eg.db)
go_non = enrichGO(eg_non$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.1, keyType = 'ENTREZID', qvalueCutoff = 0.1)

eg_neo = bitr(neo_target.nkx2_1, fromType = 'SYMBOL', toType = c('ENTREZID', 'ENSEMBL', 'SYMBOL'), OrgDb = org.Mm.eg.db)
go_neo = enrichGO(eg_neo$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.1, keyType = 'ENTREZID', qvalueCutoff = 0.1)

go_non.df = as.data.frame(go_non); go_neo.df = as.data.frame(go_neo)
go_non.df = go_non.df[order(go_non.df[,9], decreasing = T),]; go_neo.df = go_neo.df[order(go_neo.df[, 9], decreasing = T),]
## 
go_neo.df.pos = go_neo.df[grep('positive', go_neo.df$Description),]
go_neo.df.neg = go_neo.df[grep('negative', go_neo.df$Description),]


## neoplastic cells, pos
go_neo_term.pos = ggplot(go_neo.df.pos[1:8,], aes(Description, Count)) + geom_point(aes(color = pvalue, size = Count)) +
  coord_flip() + theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.y = element_text(size = 12,color = 'black'),
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        text = element_text(hjust = 0.5),
        legend.position = c(0.8, 0.4)) + scale_color_continuous(low = '#FF6D6F', high = '#4ECDC4') +
  scale_x_discrete(limits = rev(go_neo.df.pos[1:8,]$Description),
                   labels = function(x)stringr::str_wrap(x, width = 40)) +
  scale_y_continuous(breaks = seq(0, 25, 10)) +
  labs(x = 'GO Terms', y = 'Gene Numbers', title = 'Nkx2-1 targets in predicted neoplastic cells(positive regulation)')
go_neo_term.pos
ggsave(go_neo_term.pos, filename = 'Nkx2_1_target_neo_pos.pdf', width = 10, height = 16)

## neoplastic cells, neg
go_neo_term.neg = ggplot(go_neo.df.neg[1:8,], aes(Description, Count)) + geom_point(aes(color = pvalue, size = Count)) +
  coord_flip() + theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.y = element_text(size = 12,color = 'black'),
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        text = element_text(hjust = 0.5),
        legend.position = c(0.8, 0.4)) + scale_color_continuous(low = '#FF6D6F', high = '#4ECDC4') +
  scale_x_discrete(limits = rev(go_neo.df.neg[1:8,]$Description),
                   labels = function(x)stringr::str_wrap(x, width = 40)) +
  scale_y_continuous(breaks = seq(0, 25, 10)) +
  labs(x = 'GO Terms', y = 'Gene Numbers', title = 'Nkx2-1 targets in predicted neoplastic cells(negative regulation)')
go_neo_term.neg
ggsave(go_neo_term.neg, filename = 'Nkx2_1_target_neo_neg.pdf', width = 10, height = 16)
