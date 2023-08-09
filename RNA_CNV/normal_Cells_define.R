## SET UP WORKING ENVIRONMENT ##
rm(list = ls())
gc()

set.seed(1)
setwd(dir = "YOURWORKDIR")

library(dplyr)
library(scales)

source("./hclust_iterativeClustering.R")

ref_table <- read.table("./infercnv.references.txt", header = TRUE)
obs_table <- read.table("./infercnv.observations.txt", header = TRUE)
full_table <- cbind(ref_table, obs_table)

raw_expect <- MeanExpect(table = full_table,
adj_cells = colnames(ref_table), tumor_cells = colnames(obs_table))

expect_correct <- ExpectCorrection(raw_expect = raw_expect,
adj_cells = colnames(ref_table), tumor_cells = colnames(obs_table))

norm_expect_num <<- expect_correct[["norm.expect"]]
mali_expect_num <<- expect_correct[["mali.expect"]]

full_discern_K50 <- CellDiscern(table = full_table, k_group = 50, iteration = 10) # nolint
k50_pdnormal <- full_discern_K50[["pdNormal"]]
k50_pdintermediate <- full_discern_K50[["pdIntermediate"]]
k50_pdmalignant <- full_discern_K50[["pdMalignant"]]

k50 <- list(k50_pdnormal, k50_pdintermediate, k50_pdmalignant)
names(k50) <- c("pdNormal", "pdIntermediate", "pdMalignant")

saveRDS(k50, file = "./hclust_pseudobulks_k50_iterativeHclust_define.rds")
