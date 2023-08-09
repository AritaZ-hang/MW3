library(dplyr)
library(ggplot2)
#### using wilcoxon test to identify differentially expressed genes.
differentialWilcox <- function(group1, 
                              group2, 
                              dat, 
                              min.counts = 0, 
                              atac = F, 
                              diffs = 'fold change')
{
  cat("\n Filtering zero-count features across both groups...\n\n")
  group.union <- base::union(group1, group2)
  condition <- rowSums(dat[, group.union]) > min.counts
  dat <- dat[condition, ]
  print(dim(dat))
  cat(
    sprintf("Computing rank-sums on %s features ...\n",
            as.character(sum(condition)))
  )
  if(atac == T)
  {
    rns <- matrixStats::rowRanges(atac.sce)$name
    rns <- rns[condition]
  }
  else
    rns <- rownames(dat)
  out <- do.call(rbind, lapply(seq_len(nrow(dat)),
                              function(x){
                                if(x %% 10 == 0)
                                {
                                  cat(paste0(c("\t", as.character(x), " features\n")))
                                }
                                gp1 <- dat[x, group1]
                                gp2 <- dat[x, group2]
                                avg1 <- mean(gp1)
                                avg2 <- mean(gp2)
                                wt <- wilcox.test(gp1, gp2)
                                if(diffs == "fold change")
                                {
                                  diff <- log2(avg1 / avg2)
                                  out <- data.frame(
                                    gene.id = rns[x],
                                    W.statistic = wt$statistic,
                                    group.1.avg = avg1,
                                    group.2.avg = avg2,
                                    avg.log2.FC = diff,
                                    p.value = wt$p.value,
                                    row.names = F,
                                    stringsAsFactors = F
                                  )
                                }else{
                                  diff <- avg1 - avg2
                                  out <- data.frame(
                                    gene.id = rns[x],
                                    W.statistic = wt$statistic,
                                    group.1.avg = avg1,
                                    group.2.avg = avg2,
                                    differential = diff,
                                    p.value = wt$p.value,
                                    row.names = F,
                                    stringsAsFactors = F
                                  )
                                }
                              }))
  out[["p.adj"]] <- p.adjust(out[["p.value"]], method = "BH")
  return(out)
}
#### draw DEGs' Volcano plot
volcanoPlot <- function(df,
                       log.var = "avg.log2.FC",
                       pval.var = "p.adj",
                       state.var = "state",
                       do.neg.log10 = FALSE,
                       other.gene = NULL, 
                       abs_xintercept = 0.25,
                       yintercept = 0.05)
{
  df <- df[-is.infinite(df[[log.var]]), ]
  dff1 <- df %>% filter(!! sym(log.var) <= -0.5) %>% 
    slice_min(!! sym(pval.var), n = 5)
  dff2 <- df %>% filter(!! sym(log.var) >= 0.5) %>%
    slice_min(!! sym(pval.var), n = 5)
  dff <- rbind(dff1, dff2)
  other_df <- df[df$gene.id %in% other.gene,]
  dff <- rbind(dff, other_df)
  if(do.neg.log10 == TRUE)
    df[[pval.var]] <- -log10(df[[pval.var]])
  p1 <- ggplot(
    df, aes_string(x = log.var, y = pval.var)) +
    geom_point(aes_string(color = state.var), size = 1) +
    ggrepel::geom_text_repel(data = subset(df, gene.id %in% dff$gene.id), aes(label = gene.id), size = 4,
                             fontface = 'italic', box.padding = 0.2, point.padding = 0.3, segment.size = 0.3, max.overlaps = 100) + 
    scale_color_manual(values = c("#006b7b","grey", "#ef1828")) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_classic() +
    geom_vline(xintercept=c(-abs_xintercept,abs_xintercept),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(yintercept),lty=4,col="black",lwd=0.8)
  p1
  return(p1)
}
