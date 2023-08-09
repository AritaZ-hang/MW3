## SET UP WORKING ENVIRONMENT
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

## STEP 01: Prepare data ##

### load MW3 brain
MW3 = readRDS('./MW3_brain.rds')
### load unique genes file(obtained by setdiff)
unique_genes = readRDS('./unique_genes.rds')
### remove ribosomal protein genes
unique_genes = unique_genes[-grep('(Rps)|(Rpl)|(Mrps)|(Mrpl)', unique_genes)] #12571

## STEP 02: Find markers on unique genes ##

markers = FindAllMarkers(MW3, assay = 'RNA', features = unique_genes, only.pos = F)
markers.use = markers %>%
  filter(abs(avg_log2FC) > 0.5)

## STEP 03: Summarize ##

brain_types = data.frame(anno = as.character(unique(MW3$anno)))
brain_types$anno
brain_types$group = c('Glia', 'Neuron', 'Glia', 'Neuron', 'Neuron', 'Non-neuron', 'Glia', 'Glia', 'Glia', 'Non-neuron', 'Glia', 'Glia', 'Neuron', 'Glia', 'Neuron', 'Neuron', 'Neuron', 'Neuron')
names(brain_types)[1] = 'cluster'

markers = inner_join(markers, brain_types, by = 'cluster')
markers.use = inner_join(markers.use, brain_types, by = 'cluster')
table(markers.use$group)

df = data.frame(gene = markers.use$gene,
                group = markers.use$group)

df = unique(df)
table(df$group)

frac = as.data.frame(table(df$group))

## STEP 04: Summarize poly-A genes ##
### We downloaded the Poly-A genes list at VASA-seq's sTable.
vasa_polya = readRDS('./vasa_polyA_genes.rds')
table(df$gene %in% vasa_polya$gene)
names(vasa_polya)[1] = 'gene'
df_polya = inner_join(df, vasa_polya, by = 'gene')
table(df_polya$polyA, df_polya$group)
df_polya_b = as.data.frame(table(df_polya$polyA, df_polya$group))