rm(list = ls())
gc()
setwd("YOURWORKDIR")
library(tidyverse)
library(data.table)
library(stringr)
library(dplyr)
library(magrittr)
source("./hclust_pseudobulk_mapping.R")

############ 100w bin ###############
load(file = "./tumor_pseudobulks_100wbin.RData")
load(file = "./adj_pseudobulks_100wbin.RData")
load(file = "./wt_pseudobulks_100wbin.RData")

adj_pseudobulks_out.id.mapping <- generateIdMapping(adj_target_pseudobulks, type = "adj")
wt_pseudobulks_out.id.mapping <- generateIdMapping(wt_target_pseudobulks, type = "wt")
tumor_pseudobulks_out.id.mapping <- generateIdMapping(tumor_target_pseudobulks, type = "tumor")


adj_target_idx <- generateTargetIdx(adj_target_pseudobulks, adj_pseudobulks_out.id.mapping)
tumor_target_idx <- generateTargetIdx(tumor_target_pseudobulks, tumor_pseudobulks_out.id.mapping)
wt_target_idx <- generateTargetIdx(wt_target_pseudobulks, wt_pseudobulks_out.id.mapping)

saveRDS(adj_target_idx, file = "./adj_target_idx_100wbin.rds")
saveRDS(tumor_target_idx, file = "./tumor_target_idx_100wbin.rds")
saveRDS(wt_target_idx, file = "./wt_target_idx_100wbin.rds")

adj_target_idx_100wbin_merge <- do.call(rbind, adj_target_idx_100wbin)
tumor_target_idx_100wbin_merge <- do.call(rbind, tumor_target_idx_100wbin)
wt_target_idx_100wbin_merge <- do.call(rbind, wt_target_idx_100wbin)

target_idx_100wbin_merge <- rbind(adj_target_idx_100wbin_merge, tumor_target_idx_100wbin_merge, wt_target_idx_100wbin_merge)

write.csv(target_idx_100wbin_merge, file = "./target_idx_100wbin_merge.csv", row.names = F)

pseudobulks_names_100wbin <- unique(target_idx_100wbin_merge$pseudobulks_name)

write.table(pseudobulks_names_100wbin, file = "./pseudobulks_names_100wbin.tsv", sep = "\t", quote = F, row.names = F, col.names = F)
