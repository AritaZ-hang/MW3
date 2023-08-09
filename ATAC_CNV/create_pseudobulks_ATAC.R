library(data.table)
library(stringr)
library(dplyr)
set.seed(1)
setwd("YOURWORKDIR")
source("./hclust_pseudobulk_building.R")

lung_100wbin <- as.data.frame(fread(file = "./lung.wt.tumor.sort.100w.tsv"))

rownames(lung_100wbin) <- lung_100wbin$Cell_id
lung_100wbin <- lung_100wbin[, 2:dim(lung_100wbin)[2]]
lung_100wbin.t <- as.data.frame(t(as.matrix(lung_100wbin)))
metadata <- as.data.frame(fread(file = "./Lung_meta_anno.csv"))
lung_100wbin.t.use <- lung_100wbin.t

####################### pseudocell-preprocessing ################
metadata_adj <- metadata[which(metadata[["tissue"]] == "Adj"), ]
metadata_tumor <- metadata[which(metadata[["tissue"]] == "Tumor"), ]
metadata_wt <- metadata[which(metadata[["tissue"]] == "WT"), ]

rownames(metadata_adj) <- metadata_adj$cells
rownames(metadata_tumor) <- metadata_tumor$cells
rownames(metadata_wt) <- metadata_wt$real_name

dge_adj <- lung_100wbin.t.use[, metadata_adj$real_name]
dge_tumor <- lung_100wbin.t.use[, metadata_tumor$real_name]
dge_wt <- lung_100wbin.t.use[, metadata_wt$real_name]

adj_targetPseudobulks_num <- getTargetPseudobulk_num(metadata_adj, kn_wanted = 10, annotation_key = "celltype")
tumor_targetPseudobulks_num <- getTargetPseudobulk_num(metadata_tumor, kn_wanted = 100, annotation_key = "celltype")
wt_targetPseudobulks_num <- getTargetPseudobulk_num(metadata_wt, kn_wanted = 10, annotation_key = "celltype")

adj_pseudobulks_total_num <- 0
for(i in seq_along(adj_targetPseudobulks_num))
{
  adj_pseudobulks_total_num <- adj_pseudobulks_total_num + adj_targetPseudobulks_num[[i]]
}
tumor_pseudobulks_total_num <- 0
for(i in seq_along(tumor_targetPseudobulks_num))
{
  tumor_pseudobulks_total_num <- tumor_pseudobulks_total_num + tumor_targetPseudobulks_num[[i]]
}
wt_pseudobulks_total_num <- 0
for(i in seq_along(wt_targetPseudobulks_num))
{
  wt_pseudobulks_total_num <- wt_pseudobulks_total_num + wt_targetPseudobulks_num[[i]]
}

print("Adj pseudobulks total num:")
print(adj_pseudobulks_total_num)

print("tumor pseudobulks total num:")
print(tumor_pseudobulks_total_num)

print("wt pseudobulks total num:")
print(wt_pseudobulks_total_num)

print("Total pseudobulks:")
print(adj_pseudobulks_total_num + tumor_pseudobulks_total_num + wt_pseudobulks_total_num)

adj_target_pseudobulks <- generatePseudobulk(dge = dge_adj, 
                                            meta = metadata_adj, 
                                            annotation_key = "celltype", 
                                            target_num_list = adj_targetPseudobulks_num, 
                                            type = "adj")

tumor_target_pseudobulks <- generatePseudobulk(dge = dge_tumor, 
                                              meta = metadata_tumor ,
                                              annotation_key = "celltype", 
                                              target_num_list = tumor_targetPseudobulks_num,
                                              type = "tumor")

wt_target_pseudobulks <- generatePseudobulk(dge = dge_wt, 
                                           meta = metadata_wt, 
                                           annotation_key = "celltype", 
                                           target_num_list = wt_targetPseudobulks_num, 
                                           type = "wt")

total_adj <- adj_target_pseudobulks[[1]][["pseudobulk"]]
for(i in 2:length(adj_target_pseudobulks))
{
  total_adj <- cbind(total_adj, adj_target_pseudobulks[[i]][["pseudobulk"]])
}

total_tumor <- tumor_target_pseudobulks[[1]][["pseudobulk"]]
for(i in 2:length(tumor_target_pseudobulks))
{
  total_tumor <- cbind(total_tumor, tumor_target_pseudobulks[[i]][["pseudobulk"]])
}

total_wt <- wt_target_pseudobulks[[1]][['pseudobulk']]
for(i in 2:length(wt_target_pseudobulks))
{
  total_wt <- cbind(total_wt, wt_target_pseudobulks[[i]][['pseudobulk']])
}

total <- cbind(total_adj, total_tumor, total_wt)
saveRDS(total, file = "./total_pseudocell_100wbin.rds")
save(adj_target_pseudobulks, adj_targetPseudobulks_num, file = "./adj_pseudobulks_100wbin.RData")
save(tumor_target_pseudobulks, tumor_targetPseudobulks_num, file = "./tumor_pseudobulks_100wbin.RData")
save(wt_target_pseudobulks, wt_targetPseudobulks_num, file = "./wt_pseudobulks_100wbin.RData")
