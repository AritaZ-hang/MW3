## SET UP WORKING ENVIRONMENT ##
rm(list = ls())
setwd('YOURWORKDIR')
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggprism)
library(clusterProfiler)
library(org.Mm.eg.db)

## STEP 01: Load brain data from 3 kinds of technologies ## 

### load MW3
MW3 = readRDS('./MW3_brain.rds')
### load 10X
X10 = readRDS('./X10_brain.rds')
### load MCA
MCA = readRDS('./MCA_brain.rds')
### load Tabula
Tabula = readRDS('./Tabula_brain.rds')

## STEP 02: Extract count matrices and determine shared genes ##
counts_mw3 = GetAssayData(MW3, slot = 'counts', assay = 'RNA')
counts_10x = GetAssayData(X10, slot = 'counts', assay = 'RNA')
counts_mca = GetAssayData(MCA, slot = 'counts', assay = 'RNA')
counts_tabula = GetAssayData(Tabula, slot = 'counts', assay = 'RNA')

shared_genes = intersect(rownames(counts_mw3), rownames(counts_10x)) %>% intersect(rownames(counts_mca)) %>% intersect(rownames(counts_tabula))
length(shared_genes) # 14572 genes shared in 4 techniques.

## STEP 03: Identify differentially expressed genes in MW3 from other 3'-seq techniques ##
### problem is that there're no enough y observations, so a feasible way is to split the count table into several pseudo-datasets.
### We've done this by random subsampling.
### Cell-type level differential expression analysis was performed in similar procedure.

#### preprocessing...

index_1 = sample(colnames(counts_mw3), ncol(counts_mw3)*.5); counts_mw3_sample1 = counts_mw3[, index_1]
index_2 = sample(colnames(counts_mw3), ncol(counts_mw3)*.7); counts_mw3_sample2 = counts_mw3[, index_2]
index_3 = sample(colnames(counts_mw3), ncol(counts_mw3)*.3); counts_mw3_sample3 = counts_mw3[, index_3]

norm_10X = counts_10x %>% rowSums() %>% {1e6 *. / sum(.)}
norm_mca = counts_mca %>% rowSums() %>% {1e6 *. / sum(.)}
norm_tabula = counts_tabula %>% rowSums() %>% {1e6 *. / sum(.)}
norm_mw3_sample1 = counts_mw3_sample1 %>% rowSums() %>% {1e6 *. / sum(.)}
norm_mw3_sample2 = counts_mw3_sample2 %>% rowSums() %>% {1e6 *. / sum(.)}
norm_mw3_sample3 = counts_mw3_sample3 %>% rowSums() %>% {1e6 *. / sum(.)}


rna_tpm = tibble(
  "MW3_sample1" = norm_mw3_sample1[shared_genes],
  "MW3_sample2" = norm_mw3_sample2[shared_genes],
  "MW3_sample3" = norm_mw3_sample3[shared_genes],
  "10X" = norm_10X[shared_genes],
  "MCA" = norm_mca[shared_genes],
  "Tabula" = norm_tabula[shared_genes],
) %>%
  filter(MW3_sample1 != 0, `10X` != 0, MCA != 0, Tabula != 0, MW3_sample2 != 0, MW3_sample3 != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

#### two-sided T test
tpm = as.data.frame(rna_tpm); rownames(tpm) = names(rna_tpm$MW3_sample1)
types = c(rep('Total',3), rep('3-prime',3))
tpm_b = tpm
pval = c()
for(i in 1:dim(tpm)[1])
{
  p=t.test(x = tpm[i,1:3], y = tpm[i, 4:6])$p.value
  pval=c(pval,p)
}
tpm$pvalue = pval

#### calculate Log2FC
tpm$log2FC = ''
for(i in 1:dim(tpm)[1])
{
  tpm$log2FC[i] = log2(mean(sum(tpm[i, 1:3]))) - log2(mean(sum(tpm[i, 4:6])))
}

diff_genes = data.frame(log2FC = as.numeric(tpm$log2FC), 
                pvalue = as.numeric(tpm$pvalue))
rownames(diff_genes) = rownames(tpm)

#### define gene's state
diff_genes$state = ifelse(abs(diff_genes$log2FC) >= 2 & df$pvalue < 0.01, ifelse(diff_genes$log2FC < -2 & diff_genes$pvalue < 0.01, 'Down', 'Up'), 'Stable')
diff_genes$gene = rownames(diff_genes)

## STEP 04: Compute up-regulated genes' module scores in MW3 data ##

cell_types = data.frame(anno = as.character(unique(MW3$anno)))
cell_types$group = c('Glia', 'Neuron', 'Glia', 'Neuron', 'Neuron', 'Non-neuron', 'Glia', 'Glia', 'Glia', 'Non-neuron', 'Glia', 'Glia', 'Neuron', 'Glia', 'Neuron', 'Neuron', 'Neuron', 'Neuron')

brain.meta = MW3@meta.data
brain.meta = inner_join(brain.meta, cell_types, by = 'anno')
MW3 = AddMetaData(MW3, brain.meta$group, col.name = 'group')

up_genes = list(diff_genes[diff_genes$state == 'Up',]$gene)
MW3 = AddModuleScore(MW3, features = up_genes, name = 'up_genes')

## STEP 05: Annotate up-regulated genes' biotype and perform GO:BP enrichment analysis ##

gtf = as.data.frame(rtracklayer::import('mm10.gtf'))
table(up_genes[[1]] %in% gtf$gene_name) 
ups_gtf = data.frame(gene_name = up_genes[[1]])
ups_gtf$biotype = gtf[match(ups_gtf$gene_name, gtf$gene_name),]$gene_biotype
eg_ups = bitr(up_genes[[1]], fromType = 'SYMBOL', toType = c('ENTREZID'), OrgDb = 'org.Mm.eg.db')
go_ups = enrichGO(eg_ups$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP', pAdjustMethod = 'BH', qvalueCutoff = 0.05, keyType = 'ENTREZID')

