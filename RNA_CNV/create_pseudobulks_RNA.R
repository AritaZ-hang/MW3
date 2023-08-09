################## 0726 version ###################
### This time we choose kn_wanted as 200.

rm(list = ls())
setwd("YOURWORKDIR")

library(Seurat)
library(data.table)
library(multivariance)
library(fastcluster)

source("~/hclust_pseudobulk_building.R")

seob <- readRDS(file = "./seob_anno.rds")
seob_adj <- subset(seob, `tissue` == "Adj")
seob_ca <- subset(seob, `tissue` == "Tumor")
dge_adj <- GetAssayData(seob_adj@assays$RNA, slot = "counts")
dge_ca <- GetAssayData(seob_ca@assays$RNA, slot = "counts")
meta_adj <- seob_adj@meta.data
meta_ca <- seob_ca@meta.data

adj_targetPseudobulks_num <- getTargetPseudobulk_num(meta = meta_adj,
kn_wanted = 200, annotation_key = "anno")
ca_targetPseudobulks_num <- getTargetPseudobulk_num(meta = meta_ca,
kn_wanted = 200, annotation_key = "new_anno")

adj_pseudobulks_total_num <- 0
for(i in seq_along(adj_targetPseudobulks_num))
{
  adj_pseudobulks_total_num <- adj_pseudobulks_total_num + adj_targetPseudobulks_num[[i]]
}
ca_pseudobulks_total_num <- 0
for(i in seq_along(ca_targetPseudobulks_num))
{
  ca_pseudobulks_total_num <- ca_pseudobulks_total_num + ca_targetPseudobulks_num[[i]]
}

print("Adj pseudobulks total num:")
print(adj_pseudobulks_total_num)

print("CA pseudobulks total num:")
print(ca_pseudobulks_total_num)

print("Total pseudobulks:")
print(adj_pseudobulks_total_num + ca_pseudobulks_total_num)


names(adj_targetPseudobulks_num)

adj_target_pseudobulks <- generatePseudobulk(dge = dge_adj,
                                            meta = meta_adj,
                                            annotation_key = "anno",
                                            target_num_list = adj_targetPseudobulks_num,
                                            type = "adj")

ca_target_pseudobulks <- generatePseudobulk(dge = dge_ca, 
                                           meta = meta_ca ,
                                           annotation_key = "anno", 
                                           target_num_list = ca_targetPseudobulks_num,
                                           type = "ca")

total_adj <- adj_target_pseudobulks[[1]][["pseudobulk"]]
for(i in 2:length(adj_target_pseudobulks))
{
  total_adj <- cbind(total_adj, adj_target_pseudobulks[[i]][["pseudobulk"]])
}

total_ca <- ca_target_pseudobulks[[1]][["pseudobulk"]]
for(i in 2:length(ca_target_pseudobulks))
{
  total_ca <- cbind(total_ca, ca_target_pseudobulks[[i]][["pseudobulk"]])
}

total <- cbind(total_adj, total_ca)
saveRDS(total, file = ".total_pseudocell_kn200.rds")


save(adj_target_pseudobulks, adj_targetPseudobulks_num, file = "./adj_pseudobulks.RData")
save(ca_target_pseudobulks, ca_targetPseudobulks_num, file = "./ca_pseudobulks.RData")
