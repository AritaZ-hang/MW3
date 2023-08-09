rm(list = ls())
gc()
set.seed(1)
setwd(dir = "YOURWORKDIR")
library(data.table)
library(dplyr)
library(scales)
library(stringr)
library(gtools)
library(GenomicRanges)

obs_table <- as.data.frame(fread(file = "./infercnv.observations.txt"))
ref_table <- as.data.frame(fread(file = "./infercnv.references.txt"))
rownames(obs_table) <- obs_table$V1
obs_table <- obs_table[, 2:dim(obs_table)[2]]
rownames(ref_table) <- ref_table$V1
ref_table <- ref_table[, 2:dim(ref_table)[2]]

k50 <- readRDS(file = "./hclust_pseudobulks_k50_iterativeHclust_define.rds")
pdnormal <- K50$pdNormal
pdmalignant <- K50$pdMalignant
pdintermediate <- K50$pdIntermediate

processMatrix_delndup <- function(table)
{
    cnv_score_mat <- as.matrix(table)
    cnv_score_table <- as.matrix(table)
    cnv_score_table <- cnv_score_table %>% t() %>% scale() %>% rescale(to = c(-2, 2)) %>% t()
    res <- list(cnv_score_mat, cnv_score_table)
    names(res) <- c('mat', 'table')
    return(res)
}
total <- cbind(obs_table, ref_table)
res <- processMatrix_delndup(total)
pdnotnormal_table <- res$table[, c(pdintermediate, pdmalignant)]
pdnormal_table <- res$table[, pdnormal]

genes_sums_notnormal <- data.frame(sums = rowSums(pdnotnormal_table) / ncol(pdnotNormal_table))
genes_sums_normal <- data.frame(sums = rowSums(pdnormal_table)/ ncol(pdNormal_table))
genes_sums_corrected <- data.frame(sums = genes_sums_notnormal$sums - genes_sums_Normal$sums)
rownames(genes_sums_corrected) <- rownames(genes_sums_Normal)
##########################################################################
gtf = as.data.frame(rtracklayer::import("/media/ggj/Guo-4T-C2/microscopy/use/Mus_musculus.GRCm38.88.gtf"))
genes_gtf <- gtf[which(gtf$gene_name %in% rownames(genes_sums_notnormal) && gtf$type == "gene"),]
genes_gtf <- genes_gtf[!duplicated(genes_gtf$gene_name), ]
cytoband_range_ls <- readRDS(file = "./cytoband_range.rds")
############################ create overlap #############################
granges_genes <- GRanges(seqnames = genes_gtf$seqnames,
ranges = IRanges(start = genes_gtf$start, end = genes_gtf$end))
granges_cytoband <- GRanges(seqnames = cytoband_range_ls$chr,
ranges = IRanges(start = cytoband_range_ls$start, end = cytoband_range_ls$end))
overlap_result <- findOverlaps(granges_genes, granges_cytoband)
subject_hits <- subjectHits(overlap_result)
subject_hits_df <- data.frame(subject_hits = subject_hits,
                             idx = seq_along(subject_hits))
overlap_count <- countOverlaps(granges_genes, granges_cytoband, type = "any")
genes_overlap <- data.frame(genes = genes_gtf$gene_name,
                           overlap = overlap_count)
genes_overlap$idx <- cumsum(genes_overlap$overlap)
genes_overlap <- inner_join(genes_overlap, subject_hits_df, by = "idx")
############################ genes to cytoband #########################
genes2cytoband <- list()
for(i in 1:max(subject_hits))
{
  res <- list()
  cytoband <- cytoband_range_ls[["cytoband"]][i]
  genes <- genes_overlap$genes[genes_overlap$subject_hits == i]
  res[["cytoband"]] <- cytoband
  res[["genes"]] <- genes
  genes2cytoband[[i]] <- res
}
names(genes2cytoband) <- cytoband_range_ls$cytoband[1:max(subject_hits)]
############################ genes score to cytoband ####################
scores2cytoband <- list()
for(i in names(genes2cytoband))
{
  cytoband <- i
  genes <- genes2cytoband[[cytoband]][["genes"]]
  scores <- mean(na.omit(genes_sums_corrected[genes, "sums"]))
  scores2cytoband[[cytoband]] <- scores
}
scores2cytoband_df <- as.data.frame(do.call(rbind, scores2cytoband))
left_cytobands_rna <- rownames(scores2cytoband_df)[!is.nan(scores2cytoband_df$V1)]
left_cytobands_rna_scores <- scores2cytoband_df[left_cytobands_rna, ]
scores2cytoband_df <- data.frame(V1 = left_cytobands_rna_scores)
rownames(scores2cytoband_df) <- left_cytobands_rna

saveRDS(scores2cytoband_df, file = "./rna_cytoband_level_CNV_scores.rds")
