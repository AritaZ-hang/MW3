rm(list = ls())
setwd("YOURWORKDIR")
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)

source("./difgenes_wilcox_funcs.R")

seob_mw3 <- readRDS(file = "./brain_seob.rds")
seob_10x <- readRDS(file = "./10X_nuclei/brain_seob.rds")
seob_mw3_sub <- subset(seob_mw3, `celltype` %in% c("Oligodendrocytes", "Interneurons", "Purkinje neurons", "Astrocytes"))
seob_10x_sub <- subset(seob_10x, `celltype` %in% c("Oligodendrocytes", "Interneurons", "Purkinje neurons", "Astrocytes"))

counts_mw3 <- GetAssayData(seob_mw3_sub, slot = "counts")
counts_10x <- GetAssayData(seob_10x_sub, slot = "counts")

intersect_genes <- intersect(rownames(counts_mw3), rownames(counts_10x))
counts_mw3 <- as.matrix(counts_mw3[intersect_genes, ])
counts_10x <- as.matrix(counts_10x[intersect_genes, ])
counts_total <- cbind(counts_mw3, counts_10x)

meta_mw3 <- seob_mw3_sub@meta.data
meta_10x <- seob_10x_sub@meta.data
meta_mw3[["celltype"]] <- paste0(meta_mw3[["celltype"]], "_MW3")
meta_10x[["celltype"]] <- paste0(meta_10x[["celltype"]], "_10X")
meta_total <- rbind(meta_mw3, meta_10x)

merge <- CreateSeuratObject(counts = counts_total, meta.data = meta_total, min.cells = 0, min.features= 0 )
merge$celltypes <- factor(merge$celltypes)
merge@active.ident <- merge$celltypes

merge <- NormalizeData(object = merge, normalization.method = "LogNormalize", scale.factor = 10000)
merge$group <- "all"
avgr <- as.data.frame(AverageExpression(merge, group.by = "celltype")[[1]])
avgr_all <- as.data.frame(AverageExpression(merge, group.by = "group")[[1]])

retain_genes <- rownames(avgr)[which(avgr$Astrocytes_10X !=0 && avgr$Astrocytes_MW3 !=0 && avgr$Interneurons_10X !=0 && avgr$Interneurons_1MW3 !=0 && avgr$Oligodendrocytes_10X !=0 && avgr$Oligodendrocytes_MW3 !=0 && avgr$`Purkinje neurons_10X` != 0 && avgr$`Purkinje neurons_MW3` !=0)]
merge_sub <- merge[retain_genes, ]

#### sample ####
diff_astrocytes <- differentialWilcox(group1 = rownames(merge_sub@meta.data)[which(merge_sub$celltypes == "Astrocytes_MW3")], group2 = rownames(merge_sub@meta.data)[which(merge_sub$celltypes == "Astrocytes_10X")], dat = merge_sub@assays$RNA@data, min.counts = 5, atac = F, diffs = "fold_change")
diff_astrocytes$state = ifelse(abs(diff_astrocytes$avg.log2.FC) >= 2 & diff_astrocytes$p.adj <= 0.01, ifelse(diff_astrocytes$avg.log2.FC <= -2 & diff_astrocytes$p.adj <= 0.01, "Down", "Up"), "Stable")