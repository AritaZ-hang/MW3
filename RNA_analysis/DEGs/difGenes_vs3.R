rm(list = ls())
gc()
setwd("YOURWORKDIR")
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggprism)
library(clusterProfiler)
library(org.Mm.eg.db)

seob_mw3 <- readRDS("./brain_seob.rds")
seob_10x <- readRDS("./10X/brain_seob.rds")
seob_tabula <- readRDS("./tabula/brain_seob.rds")
seob_mca <- readRDS("./MCA/brain_seob.rds")

counts_mw3 <- GetAssayData(seob_mw3, slot = "counts" assay = "RNA")
counts_10x <- GetAssayData(seob_10x, slot = "counts", assay = "RNA")
counts_mca <- GetAssayData(seob_mca, slot = "counts", assay = "RNA")
counts_tabula <- GetAssayData(seob_tabula, slot = "counts", assay = "RNA")

shared_genes <- intersect(rownames(counts_mw3), rownames(counts_10x)) %>% intersect(rownames(counts_mca)) %>% intersect(rownames(counts_tabula))

## problem is that there're no enough y observations, so a feasible way is to split the count table into several pseudo-datasets.
## Also there're multiple ways: random subsampling cells; split the dataset like machine learning often does.

# sample 1: 0.5, sample 2: 0.7, sample3: 0.3

index_1 <- sample(colnames(counts_mw3), ncol(counts_mw3)*.5)
counts_mw3_sample1 <- counts_mw3[, index_1]
index_2 <- sample(colnames(counts_mw3), ncol(counts_mw3)*.7)
counts_mw3_sample2 <- counts_mw3[, index_2]
index_3 <- sample(colnames(counts_mw3), ncol(counts_mw3)*.3)
counts_mw3_sample3 <- counts_mw3[, index_3]

norm_10X <- counts_10x %>% rowSums() %>% {1e6 *. / sum(.)}
norm_mca <- counts_mca %>% rowSums() %>% {1e6 *. / sum(.)}
norm_tabula <- counts_tabula %>% rowSums() %>% {1e6 *. / sum(.)}
norm_mw3_sample1 <- counts_mw3_sample1 %>% rowSums() %>% {1e6 *. / sum(.)}
norm_mw3_sample2 <- counts_mw3_sample2 %>% rowSums() %>% {1e6 *. / sum(.)}
norm_mw3_sample3 <- counts_mw3_sample3 %>% rowSums() %>% {1e6 *. / sum(.)}

rna_tpm <- tibble(
  "MW3_sample1" = norm_mw3_sample1[shared_genes],
  "MW3_sample2" = norm_mw3_sample2[shared_genes],
  "MW3_sample3" = norm_mw3_sample3[shared_genes],
  "10X" = norm_10X[shared_genes],
  "MCA" = norm_mca[shared_genes],
  "Tabula" = norm_tabula[shared_genes],
) %>%
  filter(MW3_sample1 != 0, `10X` != 0, MCA != 0, Tabula != 0, MW3_sample2 != 0, MW3_sample3 != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

## two-sided T test & calculate log2FC.
tpm <- data.frame(
  'MW3_sample1' = rna_tpm$MW3_sample1,
  "MW3_sample2" = rna_tpm$MW3_sample2,
  "MW3_sample3" = rna_tpm$MW3_sample3,
  '10X' = rna_tpm$`10X`,
  'MCA' = rna_tpm$MCA,
  'Tabula' = rna_tpm$Tabula)
rownames(tpm) <- names(rna_tpm$MW3_sample1)

types <- c(rep('Total',3), rep('3-prime',3))

## now we can calculate the T-test p.value & log2FC.
tpm_b <- tpm
pval <- c()
for(i in 1:dim(tpm)[1])
{
  p <- t.test(x = tpm[i,1:3], y = tpm[i, 4:6])$p.value
  pval <- c(pval,p)
}
tpm$pvalue <- pval

## calculate log2FC

tpm$log2FC <- ''
for(i in 1:dim(tpm)[1])
{
  tpm$log2FC[i] <- log2(mean(sum(tpm[i, 1:3]))) - log2(mean(sum(tpm[i, 4:6])))
}

df <- data.frame(log2FC = as.numeric(tpm$log2FC), 
                pvalue = as.numeric(tpm$pvalue))
rownames(df) <- rownames(tpm)

df$state <- ifelse(abs(df$log2FC) >= 2 & df$pvalue < 0.01, ifelse(df$log2FC < -2 & df$pvalue < 0.01, 'Down', 'Up'), 'Stable')

df$gene <- rownames(df)
write.csv(df, file = './brain_mw3_3prime_diff_genes.csv', row.names = F)

library(ggrepel)
p1 = ggplot(
  df, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(aes(color = state), size = 1) +
  scale_color_manual(values = c("#006b7b","grey", "#ef1828")) +
  geom_text_repel(
    data = subset(df, pvalue < 0.001 & abs(df$log2FC) >= 3),
    aes(label = gene),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  theme(legend.position = "bottom") + 
  theme_classic() +
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 2,lty=4,col="black",lwd=0.8)

p1

ggsave(p1, filename = './total_3_prime_volcano.pdf', width = 6, height = 6)

## group these into neurons, glias & non-neurons.
brain_types <- data.frame(celltypes = as.character(unique(seob_mw3$celltypes)))
brain_types$group <- c('Neuron', "Neuron", "Non-neuron", "Neuron", "Glia", "Glia", "Neuron", "Glia", "Non-neuron", "Glia", "Neuron", "Neuron", "Non-neuron", "Neuron")

brain.meta <- seob_mw3@meta.data
brain.meta <- inner_join(brain.meta, brain_types, by = "celltypes")
seob_mw3 <- AddMetaData(seob_mw3, brain.meta$group, col.name = "group")
up_genes <- list(df[df$state == "Up",]$gene)
seob_mw3 <- AddModuleScore(seob_mw3, features = up_genes, name = "up_genes")
down_genes <- list(df[df$state == "Down",]$gene)
seob_mw3 <- AddModuleScore(seob_mw3, features = down_genes, name = 'down_genes')

ups_df <- seob_mw3@meta.data[, c("up_genes1", "group")]
ups_df$group <- factor(ups_df$group, levels = c("Glia", "Non-neuron", "Neuron"))

colors <- brewer.pal(3, 'Set2')

p2 <- ggplot(ups_df, aes(x = group, y = up_genes1)) + 
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(0.3), aes(fill = group), color = "black") +
  geom_boxplot(width = 0.2, aes(fill = group), notch = F, outlier.shape = 0, outlier.size = 0 ,outlier.stroke = 0, position = position_dodge(0.3), outlier.alpha = 0) +
  theme_classic() + 
  geom_signif(comparisons = list(c("Glia", "Non-neuron"), c("Non-neuron", "Neuron")), map_signif_level = F, y_position = 0.35) +
  labs(x = "Categories of Cell Types", y = "Module Scores") +
  scale_fill_manual(values = colors) + 
  theme(legend.title = element_blank())
p2
ggsave(p2, filename = "./ups_genes_mw3_scores.pdf", width = 6, height = 6)

gtf <- as.data.frame(rtracklayer::import("./Mus_musculus.GRCm38.88.gtf"))
ups_gtf <- data.frame(gene_name = up_genes[[1]])
ups_gtf$biotype <- gtf[match(ups_gtf$gene_name, gtf$gene_name),]$gene_biotype

## GO enrichment
eg_ups <- bitr(up_genes[[1]], fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Mm.eg.db")
go_ups <- enrichGO(eg_ups$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, keyType = "ENTREZID")

go_ups.use <- go_ups[c(1:20),]
go_ups.use <- go_ups.use[order(go_ups.use[,9],decreasing = TRUE),]

p3 = ggplot(go_ups.use,aes(Description,Count))+geom_point(aes(color= pvalue,size=Count))+
  labs(x = "GO:BP Terms",y = "Gene Numbers",title = "") + 
  coord_flip()+theme_bw() +
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.8, size = 15), 
        axis.text.y = element_text(size = 15,color='black'),
        axis.text.x = element_text(size = 15,color='black'), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        text = element_text(hjust = 0.5))+
        scale_color_continuous(low='#FF6D6F',high='#4ecdc4')+
  scale_x_discrete(limits = rev(go_ups.use$Description),
                   labels = function(x) stringr::str_wrap(x, width = 40))+
  theme(text = element_text(size=40))+
  scale_y_continuous(breaks=seq(0,25,10))+theme(legend.position=c(0.8,0.4))
p3

ggsave(p3, filename = './up_genes_go_bp.pdf', width = 10, height = 16)
