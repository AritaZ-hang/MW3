rm(list = ls())
gc()
library(tidyverse)
library(data.table)
library(stringr)
library(dplyr)
library(magrittr)
library(scales)
library(stringr)
library(gtools)
library(GenomicRanges)
library(ggplot2)

setwd("YOURWORKDIR")
load(file = "./Alleloscope/use_mat.RData")

#################### deal with the mat #####################

smooth_mat_norm <- ((smooth_mat_norm) %>% t() %>% scale() %>% rescale(to = c(-2, 2)) %>% t())
nmf_clusters <- as.data.frame(fread(file = "./copyscat/100wbin/lung_cnv_total_1e+06_nmf_clusters.csv"))
nmf_clusters$group <- "malignant"
nmf_clusters$group[which(nmf_clusters$nmf_results.cellAssigns == 1)] <- "normal"
nmf_clusters$Barcode <- gsub("-1", "", nmf_clusters$Barcode)

norm_mat_normal <- smooth_mat_norm[, nmf_clusters$Barcode[which(nmf_clusters$group == "normal")]]
norm_mat_malignant <- smooth_mat_norm[, nmf_clusters$Barcode[which(nmf_clusters$group == "malignant")]]

norm_mat_normal_sums <- data.frame(sums = rowSums(norm_mat_normal) / ncol(norm_mat_normal))
norm_mat_malignant_sums <- data.frame(sums = rowSums(norm_mat_malignant) / ncol(norm_mat_malignant))

norm_mat_correct <- data.frame(sums = norm_mat_malignant_sums$sums - norm_mat_normal_sums$sums)
rownames(norm_mat_correct) <- rownames(norm_mat_normal_sums)

chr_ranges <- rownames(smooth_mat_norm)
chr_ranges <- reshape2::colsplit(chr_ranges, pattern = "-", names = c("chr", "start", "end"))

cytoband_range <- readRDS(file = "./cytoband_range.rds")
chr_ranges$chr <- gsub("chr", "", chr_ranges$chr)
chr_ranges$bin <- paste("chr",chr_ranges$chr, "-", chr_ranges$start, "-", chr_ranges$end, sep = "")

#################### create overlap ######################
granges_chr <- GRanges(seqnames = chr_ranges$chr,
ranges = IRanges(start = chr_ranges$start, end = chr_ranges$end))
granges_cytoband <- GRanges(seqnames = cytoband_range$chr,
ranges = IRanges(start = cytoband_range$start, cytoband_range$end))
overlap_result <- findOverlaps(granges_chr, granges_cytoband, type = "any")
query_hits <- queryHits(overlap_result)
subject_hits <- subjectHits(overlap_result)
subject_hits_df <- data.frame(subject_hits = subject_hits,
                             idx = seq_along(subject_hits))
overlap_count <- countOverlaps(granges_chr, granges_cytoband, type = "any")
bins_overlap <- data.frame(bins = chr_ranges$bin, 
                          overlap = overlap_count)
bins_overlap$idx <- cumsum(bins_overlap$overlap)
bins_overlap <- inner_join(bins_overlap, subject_hits_df, by = "idx")
#################### bin to cytoband ######################
bins2cytoband <- list()
for(i in 1:max(subject_hits))
{
  res <- list()
  cytoband <- cytoband_range[["cytoband"]][i]
  bins <- bins_overlap$bins[bins_overlap$subject_hits == i]
  res[["cytoband"]] <- cytoband
  res[["bins"]] <- bins
  bins2cytoband[[i]] <- res
}
names(bins2cytoband) <- cytoband_range$cytoband[1:max(subject_hits)]
################## scores to cytoband ######################
alleloscope_cytoband_cnv_scores <- list()
for(i in names(bins2cytoband))
{
  cytoband <- i
  bins <- bins2cytoband[[cytoband]][["bins"]]
  scores <- mean(na.omit(norm_mat_correct[bins, "sums"]))
  alleloscope_cytoband_cnv_scores[[cytoband]] <- scores
}
alleloscope_cytoband_cnv_scores_df <- as.data.frame(do.call(rbind, alleloscope_cytoband_cnv_scores))


########### copyscat ##########
copyscat_cytoband_cnv_scores <- as.data.frame(fread(file = "./copyscat/100wbin/lung_cnv_total_1e+06_clean_cnv_cnv_scores.csv"))
rownames(copyscat_cytoband_cnv_scores) <- copyscat_cytoband_cnv_scores$rowname
copyscat_cytoband_cnv_scores <- copyscat_cytoband_cnv_scores[, 2:dim(copyscat_cytoband_cnv_scores)[2]]
copyscat_cytoband_cnv_scores <- t(copyscat_cytoband_cnv_scores %>% t() %>% scale() %>% rescale(to = c(-2, 2)) %>% t())
nmf_clusters <- as.data.frame(fread(file = "./copyscat/100wbin/lung_cnv_total_1e+06_nmf_clusters.csv"))
nmf_clusters$group <- "malignant"
nmf_clusters$group[which(nmf_clusters$nmf_results.cellAssigns == 1)] <- "normal"
copyscat_cytoband_cnv_scores_normal <- copyscat_cytoband_cnv_scores[, nmf_clusters$Barcode[which(nmf_clusters$group == "normal")]]
copyscat_cytoband_cnv_scores_notNormal <- copyscat_cytoband_cnv_scores[, nmf_clusters$Barcode[which(nmf_clusters$group == "malignant")]]
copyscat_cytoband_cnv_scores_normal_sums <- data.frame(sums = rowSums(copyscat_cytoband_cnv_scores_normal) / ncol(copyscat_cytoband_cnv_scores_normal))
copyscat_cytoband_cnv_scores_notNormal_sums <- data.frame(sums = rowSums(copyscat_cytoband_cnv_scores_notNormal) / ncol(copyscat_cytoband_cnv_scores_notNormal))
copyscat_cytoband_cnv_scores_correct_sums <- data.frame(sums = copyscat_cytoband_cnv_scores_notNormal_sums$sums - copyscat_cytoband_cnv_scores_normal_sums$sums)
rownames(copyscat_cytoband_cnv_scores_correct_sums) <- rownames(copyscat_cytoband_cnv_scores_normal)

##################################################################
intersect_cytoband <- intersect(rownames(alleloscope_cytoband_cnv_scores_df), rownames(copyscat_cytoband_cnv_score_correct_sums))
alleloscope_cytoband_cnv_scores_df_sub <- data.frame(scores = alleloscope_cytoband_cnv_scores_df[intersect_cytoband,], 
                                         cytoband = intersect_cytoband)
copyscat_cytoband_cnv_scores_correct_sums_sub <- data.frame(scores = copyscat_cytoband_cnv_scores_correct_sums[intersect_cytoband,], 
                                                 cytoband = intersect_cytoband)
ATACs <- data.frame(copyscat = copyscat_cytoband_cnv_scores_correct_sums_sub[['scores']], 
alleloscope = alleloscope_cytoband_cnv_scores_df_sub[['scores']], 
cytoband = intersect_cytoband)

RNAs <- readRDS(file = "./rna_cytoband_level_CNV_scores.rds")
names(RNAs) <- c("inferCNV", "cytoband")
total <- inner_join(ATACs, RNAs, by = "cytoband")
total$group.new <- "not correlated"
total$group.new[which(total$infercnv >=0 && total$copyscat >=0 && total$alleloscope >=0)] <- "dup effects"
total$group.new[which(total$infercnv <=0 & total$copyscat <=0 & total$alleloscope <=0)] <- "del effects"

saveRDS(total, file = "./compare_RNA_ATAC_prediction.rds")