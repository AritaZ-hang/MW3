## SET UP WORKING ENVIRONMENT ##
setwd('YOURWORKDIR')
rm(list = ls()); gc()
set.seed(1)
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

## STEP 01: Read in the enrichment score matrix produced at 'SCRIP enrich' and processing ##

enri = fread('SCRIP_enrichment.txt', header = T)
enri = as.data.frame(t(enri))
enri[1:5, 1:5]
colnames(enri) = enri[1,]
enri_b = enri[2:dim(enri)[1],]

seurat = CreateSeuratObject(counts = enri_b, project = 'brain')
seurat@assays$RNA@scale.data = as.matrix(seurat@assays$RNA@counts)
seurat = FindVariableFeatures(seurat, selection.method = 'vst', nfeatures = 2000)
seurat = RunPCA(seurat, features = VariableFeatures(seurat))

seurat = FindNeighbors(seurat, dims = 1:25)
seurat = FindClusters(seurat, resolution = 0.8)

seurat = RunUMAP(seurat, dims = 1:30)
DimPlot(seurat, reduction = 'umap')


## STEP 02: Change Seurat object 'seurat''s active.ident into your determined cell type annotations ##
### This step is skipped.

## STEP 03: Find differentially expressed TFs across cell types ##
markers = FindAllMarkers(seurat, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1)
markers.use = brain.markers %>% 
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)
gene = unique(markers.use$gene)

## STEP 04: Draw Heatmap for differentially expressed TFs ##

mat = GetAssayData(seurat, slot = 'counts')
mat_use = as.matrix(mat[gene,])
mycol = brewer.pal(length(seurat$celltype), 'Spectral') 
names = as.data.frame(seurat_b@active.ident)
colnames(names) = 'cell_type'
type = names$cell_type
ha = HeatmapAnnotation(type = type, annotation_name_side = 'left',
                       col = list(type = c("CELLTYPE1" = mycol[1], 'CELLTYPE2' = mycol[2])))

## Perform min-max scaling on the matrix
for(i in 1:dim(mat_use)[1])
{
  min = min(mat_use[i,])
  max = max(mat_use[i,])
  for(c in 1:dim(mat_use)[2])
  {
    mat_use[i,c] = (mat_use[i,c] - min) / (max-min)
  }
}
names$cell_type = as.factor(names$cell_type)
annotation_col = names

options(repr.plot.width = 10, repr.plot.height = 20, repr.plot.res = 300)
p = Heatmap(mat_use, cluster_rows = T,
            cluster_columns = T,
            show_column_names = F,
            show_row_names = T,
            top_annotation = ha,
            column_km = 12,
            row_km = 15,
            show_column_dend = F,
            show_row_dend = F,
            heatmap_legend_param = list(
              title = 'score',
              title_position = 'leftcenter-rot'
            ),
            row_gap = unit(0, 'mm'),
            column_gap = unit(0, 'mm'),
            border = T,
            width = unit(40, 'cm'),
            height = unit(75, 'cm'),
            row_names_side = 'left')
pdf(file = 'heatmap.pdf', width = 20, height = 30)
p
dev.off()
