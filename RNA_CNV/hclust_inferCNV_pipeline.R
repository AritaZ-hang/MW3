setwd("YOURWORKDIR")
rm(list = ls())
gc()
library(data.table)
library(dplyr)
library(BiocParallel)
library(infercnv)
library(scales)

pseudo <- readRDS(file = "./total_pseudocell_kn200.rds")
meta <- data.frame(barcode = colnames(pseudo))
meta$celltype <- gsub("_Pseudobulk_[0-9]+", "", meta$barcode)
meta$group <- "CA"
meta$group[grep("adj", meta$celltype)] = "Adj"

k50 <- readRDS(file = "./hclust_pseudobulks_k50_iterativeHclust_define.rds")
predicted <- data.frame(barcode = c(k50$pdNormal, k50$pdIntermediate, k50$pdMalignant), 
                        predicted = c(rep("pdNormal", length(k50$pdNormal),) rep("pdIntermediate", length(k50$pdIntermediate)), rep("pdMalignant", length(k50$pdMalignant))))
meta <- inner_join(meta, predicted, by = "barcode")

geneorder <- readRDS(file = "./gencode.vM31.pos.rds")
geneorder <- geneorder[!duplicatred(geneorder$gene_name),]
rownames(geneorder) <- generoder$gene_name

genecommon <- intersect(rownames(geneorder), rownames(pseudo))
pseudo <- pseudo[genecommon, ]
geneorder <- geneorder[genecommon,]
geneorder <- geneorder[, c("seqnames", "start", "end")]
names(geneorder) = c("chr", "start", "stop")

counts <- as.matrix(pseudo)
metause.tmp <- meta[, c("predicted")]
metause.tmp <- as.data.frame(metause.tmp)
rownames(metause.tmp) <- meta$barcode
colnames(metause.tmp) <- "predicted"

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts, 
annotations_file = metause.tmp, gene_order_file = geneorder, 
ref_group_names = "pdNormal")
options(scipen = 100)
options(expressions = 500000)
infercnv_obj <- run(infercnv_obj = infercnv_obj, 
cutoff = 1, out_dir = "./hclust_results/", cluster_by_groups = T, denoise = T, HMM = T, analysis_mode = "sample", num_threads = 16, num_ref_groups = NULL, output_format = "pdf", scale_data = F, write_expr_matrix = T)
save(infercnv_obj, file = "./hclust_results/infer_cnv_obj.RData")