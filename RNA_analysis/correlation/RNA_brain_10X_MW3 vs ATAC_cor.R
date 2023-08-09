rm(list = ls())
gc()
setwd("YOURWORKDIR")
library(Seurat)
library(ArchR)
library(reshape2)
addArchRThreads(threads = 12)
addArchRGenome("mm10")
source("./ArchR_Wrapped_funcs.R")

proj <- readRDS("./brain_ArchR.rds")
seob_mw3 <- readRDS("./brain_seob.rds")

## load 10x
seob_10x <- readRDS("./10X/brain_seob.rds")

## integration with mw3.
proj_labeltrans_mw3 <- LabelTransfer(proj, seob_mw3, anno_name = 'celltypes', res = 0.5, TileorPeak = 'Tile')
predictScore_mw3 <- proj_labeltrans_mw3$predictedScore_Un
predictScore_mw3_df <- data.frame(cells = proj_labeltrans_mw3$cellNames, 
                                 predictedScore = predictScore_mw3)
mw3_predicted_score_p <- ggplot(predictScore_mw3_df, aes(x = predictedScore)) + 
  geom_histogram(aes(y = ..density..),color = "royalblue", fill = "lightblue") + 
  geom_density(alpha = .2, fill = "darkblue") + 
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") + ggtitle("Label Transfer by MW3-RNA Brain") + theme_bw()
ggsave(mw3_predicted_score_p, filename = "./MW3_predicted_score_hist.pdf", width = 6, height = 4)

proj_labeltrans_mw3.show <- proj_labeltrans_mw3[proj_labeltrans_mw3$predictedScore_Un >= 0.5,]

prediction_df.mw3 <- data.frame(cells = proj_labeltrans_mw3.show$cellNames,
                               UMAP_1 = proj_labeltrans_mw3.show@embeddings@listData$UMAP_bin_res0.5@listData$df$`IterativeLSI_bin_res0.5#UMAP_Dimension_1`,
                               UMAP_2 = proj_labeltrans_mw3.show@embeddings@listData$UMAP_bin_res0.5@listData$df$`IterativeLSI_bin_res0.5#UMAP_Dimension_2`,
                               anno = proj_labeltrans_mw3.show$predictedGroup_Un)

pdf(file = "./labeltransfer_MW3_0.5cutoff.pdf", width = 10, height = 10)
ggPoint(x = prediction_df.mw3$UMAP_1,
        y = prediction_df.mw3$UMAP_2,
        color = prediction_df.mw3$anno,
        size = 0.3,
        labelSize  = 0,
        discreteSet = "kelly",
        colorTitle = "Predicted Cell Types") +
  theme_classic() + 
  theme(legend.position = "right",
  axis.ticks = element_blank(),
  axis.text = element_blank()) + 
  xlab("UMAP dimension 1") + 
  ylab("UMAP dimension 2")
dev.off()

save(proj_labeltrans_mw3, proj_labeltrans_mw3.show, prediction_df.mw3, file = './labeltrans_mw3_data.RData')

## integration with 10x.
proj_labeltrans_10x <- LabelTransfer(proj, seob_10X, anno_name = "anno", res = 0.5, TileorPeak = "Tile")

predictScore_10x <- proj_labeltrans_10x$predictedScore_Un
predictScore_10x_df <- data.frame(cells = proj_labeltrans_10x$cellNames, 
                                 predictedScore = predictScore_10x)
x10_predicted_score_p <- ggplot(predictScore_10x_df, aes(x = predictedScore)) + 
  geom_histogram(aes(y = ..density..),color = "royalblue", fill = "lightblue") + 
  geom_density(alpha = .2, fill = "darkblue") + 
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") + ggtitle("Label Transfer by 10X-RNA Brain") + theme_bw()
ggsave(x10_predicted_score_p, filename = "./10X_predicted_score_hist.pdf", width = 6, height = 4)

proj_labeltrans_10x.show <- proj_labeltrans_10x[proj_labeltrans_10x$predictedScore_Un >= 0.5,]
prediction_df.10x <- data.frame(cells = proj_labeltrans_10x.show$cellNames,
                               UMAP_1 = proj_labeltrans_10x.show@embeddings@listData$UMAP_bin_res0.5@listData$df$`IterativeLSI_bin_res0.5#UMAP_Dimension_1`,
                               UMAP_2 = proj_labeltrans_10x.show@embeddings@listData$UMAP_bin_res0.5@listData$df$`IterativeLSI_bin_res0.5#UMAP_Dimension_2`,
                               anno = proj_labeltrans_10x.show$predictedGroup_Un)

pdf(file = "./labeltransfer_10x_0.5cutoff.pdf", width = 10, height = 10)
ggPoint(x = prediction_df.10x$UMAP_1,
        y = prediction_df.10x$UMAP_2,
        color = prediction_df.10x$anno,
        size = 0.3,
        labelSize  = 0,
        discreteSet = "kelly",
        colorTitle = "Predicted Cell Types") + 
  theme_classic() + theme(legend.position = "right", axis.ticks = element_blank(), axis.text = element_blank()) + xlab("UMAP dimension 1") + ylab("UMAP dimension 2")
dev.off()

save(proj_labeltrans_10x, proj_labeltrans_10x.show, prediction_df.10x, file = './labeltrans_10x_data.RData')

## spearman correlation
Expr_perlineage <- function(expM, proj, anno)
{
  if(class(proj)[1] == "ArchRProject")
  {
    annos <- ArchR::getCellColData(proj, select = anno)
    annos <- as.data.frame(annos)
  }
  else if(class(proj)[1] == "Seurat")
  {
    annos <- proj@meta.data[, anno]
    annos <- as.data.frame(annos)
    rownames(annos) <- colnames(proj)
  }
  Expr_perlineage <- c()
  for( i in unique(annos[,1]))
  {
    selected.cells <- rownames(annos)[annos[,1] == i]
    expM.extracted <- expM[,selected.cells]
    avg_expr <- apply(expM.extracted, 1, mean)
    Expr_perlineage <- cbind(Expr_perlineage, avg_expr)
  }
  Expr_perlineage <- as.data.frame(Expr_perlineage)
  colnames(Expr_perlineage) <- unique(annos[,1])
  return(Expr_perlineage)
}
GetCorrelationMatrix <- function(Rna_lineage_expr, Atac_lineage_expr)
{
  row_num <- dim(Rna_lineage_expr)[2]
  col_num <- dim(Atac_lineage_expr)[2]
  CorrelationMatrix <- matrix(nrow = row_num, ncol = col_num)
  for(i in 1:row_num)
  {
    for(j in 1:col_num)
    {
      CorrelationMatrix[i,j] <- cor(Rna_lineage_expr[,i], Atac_lineage_expr[,j], method = "spearman")
    }
  }
  CorrelationMatrix <- as.data.frame(CorrelationMatrix)
  names(CorrelationMatrix) <- names(Atac_lineage_expr)
  rownames(CorrelationMatrix) <- names(Rna_lineage_expr)
  return(CorrelationMatrix)
}

HGs.mw3 <- VariableFeatures(seob_mw3)
GeneScoreMat <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
genemat <- GeneScoreMat@assays@data@listData$GeneScoreMatrix
genes <- as.character(GeneScoreMat@elementMetadata@listData[["name"]])
rownames(genemat) <- genes
HGs.mw3 <- HGs.mw3[HGs.mw3 %in% genes]
rna_expM.mw3 <- seob_10X@assays$RNA@data[HGs.mw3,]
atac_expM.mw3 <- as.data.frame(genemat[HGs.mw3,])
save(rna_expM.mw3, atac_expM.mw3, file = "./mw3_rna_atac_expM.RData")

Rna_lineage_expr.mw3 <- Expr_perlineage(rna_expM.mw3, pbmc.mw, anno = 'celltypes')
Atac_lineage_expr.mw3 <- Expr_perlineage(atac_expM.mw3, proj, anno = 'final_anno')
correlation_matrix.mw3 <- GetCorrelationMatrix(Rna_lineage_expr.mw3, Atac_lineage_expr.mw3)
write.csv(correlation_matrix.mw3, file = "./mw3_spearmanCor_matrix.csv")

library(viridis)
cor_mw3 <- pheatmap::pheatmap(correlation_matrix.mw3,
                             color = viridis(10),
                             cluster_rows = T,
                             cluster_cols = T,
                             angle_col = 315,
                             fontsize = 10,
                             border_color = "white",
                             scale = "none",
                             width = 5,
                             height = 5)

correlation_matrix.mw3 <- correlation_matrix.mw3[cor_mw3$tree_row$order, cor_mw3$tree_col$order]
pheatmap::pheatmap(correlation_matrix.mw3,
                   color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(50),
                   cluster_rows = F,
                   cluster_cols = F,
                   angle_col = 315,
                   fontsize = 10,
                   border_color = "white",
                   scale = "none",
                   width = 10,
                   height = 10,
                   filename = "./mw3_spearmanCor.pdf")
## 10x.
HGs.10x <- VariableFeatures(seob_10x)
GeneScoreMat <- getMatrixFromProject(proj, useMatrix = 'GeneScoreMatrix')
genemat <- GeneScoreMat@assays@data@listData$GeneScoreMatrix
genes <- as.character(GeneScoreMat@elementMetadata@listData[["name"]])
rownames(genemat) <- genes
HGs.10x <- HGs.10x[HGs.10x %in% genes]
rna_expM.10x <- seob_10x@assays$RNA@data[HGs.10x,]
atac_expM.10x <- as.data.frame(genemat[HGs.10x,])
save(rna_expM.10x, atac_expM.10x, file = './10x_rna_atac_expM.RData')


Rna_lineage_expr.10x <- Expr_perlineage(rna_expM.10x, seob_10x, anno = "anno")
Atac_lineage_expr.10x <- Expr_perlineage(atac_expM.10x, proj, anno = 'final_anno')


correlation_matrix.10x = GetCorrelationMatrix(Rna_lineage_expr.10x, Atac_lineage_expr.10x)
write.csv(correlation_matrix.10x, file = '/media/ggj/TOSHIBA EXT/MW3_benchmark/RNA_brain/RNA_brain_10X_MW3 vs ATAC_cor/10x_spearmanCor_matrix.csv')

library(viridis)
cor_10x <- pheatmap::pheatmap(correlation_matrix.10x,
                             color = viridis(10),
                             cluster_rows = T,
                             cluster_cols = T,
                             angle_col = 315,
                             fontsize = 10,
                             border_color = "white",
                             scale = "none",
                             width = 5,
                             height = 5)

correlation_matrix.10x <- correlation_matrix.10x[cor_10x$tree_row$order, cor_10x$tree_col$order]
pheatmap::pheatmap(correlation_matrix.10x,
                   color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(50),
                   cluster_rows = F,
                   cluster_cols = F,
                   angle_col = 315,
                   fontsize = 10,
                   border_color = "white",
                   scale = "none",
                   width = 10,
                   height = 10,
                   filename = "./10x_spearmanCor.pdf")