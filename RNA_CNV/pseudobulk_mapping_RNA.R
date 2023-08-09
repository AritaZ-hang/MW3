rm(list = ls())
setwd("YOURWORKDIR")

library(Seurat)
library(data.table)
library(dplyr)

source("./hclust_pseudobulk_mapping.R")

seob <- readRDS(file = "./seob_anno.rds")
load(file = "./adj_pseudobulks.RData")
load(file = "./ca_pseudobulks.RData")
k50 <- readRDS(file = "./hclust_pseudobulks_k50_iterativeHclust_define.rds")


adj_pseudobulks_out.id.mapping <- generateIdMapping(adj_target_pseudobulks,
type = "adj")
ca_pseudobulks_out.id.mapping <- generateIdMapping(ca_target_pseudobulks,
type = "ca")


adj_target_idx <- generateTargetIdx(adj_target_pseudobulks,
adj_pseudobulks_out.id.mapping)
ca_target_idx <- generateTargetIdx(ca_target_pseudobulks,
ca_pseudobulks_out.id.mapping)


k50_predict <- data.frame(pseudobulks_name = c(k50$pdNormal, k50$pdMalignant, k50$pdIntermediate), 
                         predict_new = c(rep("pdNormal", length(k50$pdNormal)), rep("pdMalignant", length(k50$pdMalignant)), rep("pdIntermediate", length(k50$pdIntermediate))))

total_adj <- do.call(rbind, adj_target_idx)
total_ca <- do.call(rbind, ca_target_idx)

total <- rbind(total_adj, total_ca)
total <- inner_join(total, k50_predict, by = "pseudobulks_name")

old_meta <- seob@meta.data
old_meta$cells <- rownames(old_meta)
old_meta <- left_join(old_meta, total, by = "cells")

rownames(old_meta) <- old_meta[["cells"]]
old_meta[["predict"]][is.na(old_meta[["predict_new"]])] <- "Filtered"

seob$predict_new <- old_meta$predict_new
seob$pseudobulks_name <- old_meta$pseudobulks_name

saveRDS(seob, file = "./hclust_pseudobulk_anno_seob.rds")
saveRDS(old_meta, file = "./old_meta.rds")