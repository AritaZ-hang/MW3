## construct pseudo samples.
## count_mw3[, cells] as input.
sampling = function(counts)
{
  index_1 = sample(colnames(counts), ncol(counts)*.3)
  counts_sample1 = counts[, index_1]

  index_2 = sample(colnames(counts), ncol(counts)*.5)
  counts_sample2 = counts[, index_2]

  index_3 = sample(colnames(counts), ncol(counts)*.7)
  counts_sample3 = counts[, index_3]

  list = list()
  list[[1]] = counts_sample1
  list[[2]] = counts_sample2
  list[[3]] = counts_sample3

  return(list)
}

## normalize to tpm.
tpm_norm = function(sample_list_a, sample_list_b, shared_genes)
{
  norm_a = sample_list_a[[1]] %>% rowSums() %>% {1e6 *. / sum(.)}
  norm_b = sample_list_a[[2]] %>% rowSums() %>% {1e6 *. / sum(.)}
  norm_c = sample_list_a[[3]] %>% rowSums() %>% {1e6 *. / sum(.)}
  norm_d = sample_list_b[[1]] %>% rowSums() %>% {1e6 *. / sum(.)}
  norm_e = sample_list_b[[2]] %>% rowSums() %>% {1e6 *. / sum(.)}
  norm_f = sample_list_b[[3]] %>% rowSums() %>% {1e6 *. / sum(.)}

  rna_tpm = tibble(
    "sample1_a" = norm_a[shared_genes],
    "sample2_a" = norm_b[shared_genes],
    "sample3_a" = norm_c[shared_genes],
    "sample1_b" = norm_d[shared_genes],
    "sample2_b" = norm_e[shared_genes],
    "sample3_b" = norm_f[shared_genes],
  ) %>%
    filter(sample1_a != 0, sample2_a != 0, sample3_a != 0, sample1_b != 0, sample2_b != 0, sample3_b != 0) %>%
    mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

  tpm = data.frame(
    "sample1_a" = rna_tpm$sample1_a,
    "sample2_a" = rna_tpm$sample2_a,
    "sample3_a" = rna_tpm$sample3_a,
    'sample1_b' = rna_tpm$sample1_b,
    'sample2_b' = rna_tpm$sample2_b,
    'sample3_b' = rna_tpm$sample3_b)
  rownames(tpm) = names(rna_tpm$sample1_a)

  return(tpm)
}

## get diff genes by t.test & log2FC.
diffGenes = function(tpm, p.val, log2th, type)
{
  pval = c()
  for(i in 1:dim(tpm)[1])
  {
    p=t.test(x = tpm[i,1:3], y = tpm[i, 4:6])$p.value
    pval=c(pval,p)
  }
  tpm$pvalue = pval

  for(i in 1:dim(tpm)[1])
  {
    tpm$log2FC[i] = log2(mean(sum(tpm[i, 1:3]))) - log2(mean(sum(tpm[i, 4:6])))
  }

  df = data.frame(log2FC = as.numeric(tpm$log2FC),
                  pvalue = as.numeric(tpm$pvalue))
  rownames(df) = rownames(tpm)

  df$state = ifelse(abs(df$log2FC) >= log2th & df$pvalue < p.val, ifelse(df$log2FC < -log2th & df$pvalue < p.val, 'Down', 'Up'), 'Stable')
  df$gene = rownames(df)
  df$group = type
  return(df)
}

volcano_p = function(dif_genes, polar = T, topGeneN = 5, orderby = c('log2FC'))
{
  # get background cols.
  purrr::map_df(unique(dif_genes$group), function(x)
  {
    tmp = dif_genes %>%
      dplyr::filter(group == x)
    new.tmp = data.frame(group = x,
                         min = min(tmp$log2FC) - 0.2,
                         max = max(tmp$log2FC) + 0.2)
    return(new.tmp)
  }) -> back.data

  ## get top genes
  top.marker.tmp = dif_genes %>%
    dplyr::group_by(state, group)

  ## order
  if(length(orderby == 1))
  {
    top.marker.max = top.marker.tmp %>%
      dplyr::slice_max(n = topGeneN, order_by = get(orderby))
    top.marker.min = top.marker.tmp %>%
      dplyr::group_by(group) %>%
      slice_min(n = topGeneN, order_by = get(orderby))
  }
  else
  {
    top.marker.max = top.marker.tmp %>%
      arrange(desc(get(orderby[1])), get(orderby[2])) %>%
      slice_head(n = topGeneN)
    top.marker.min = top.marker.tmp %>%
      arrange(desc(get(orderby[1])), get(orderby[2])) %>%
      slice_tail(n = topGeneN)
  }
  top.marker = rbind(top.marker.max, top.marker.min)

  ## plot
  back.col = 'grey93'
  aesCol = c('#0099CC', '#CC3333')
  p1 = ggplot(dif_genes,
              aes(x = group, y = log2FC)) +
    geom_col(data = back.data,
             aes(x = group, y = min), fill = back.col) +
    geom_col(data = back.data,
             aes(x = group, y = min), fill = back.col) +
    geom_jitter(aes(color = state), size = 0.75) +
    scale_color_manual(values = c(aesCol[2], aesCol[1]))

  ## add themes
  p2 = p1 +
    scale_y_continuous(n.breaks = 6) +
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank(),
          legend.position = 'right',
          legend.title = element_blank(),
          legend.background = element_blank()) +
    xlab('Cell types') + ylab('Log2 Fold Change') +
    guides(color = guide_legend(override.aes = list(size = 5)))

  # add tile
  p3 = p2 +
    geom_tile(aes(x = group, y = 0, fill = group),
              color = 'black',
              height = 2*2,
              alpha = 0.3,
              show.legend = F) +
    ggrepel::geom_label_repel(data = top.marker,
                              aes(x = group, y = log2FC, label = gene),
                              max.overlaps = 50)
  ## polar
  if(polar == TRUE)
  {
    p4 = p3 +
      scale_y_continuous(n.breaks = 6,
                         expand = expansion(mult = c(-1,1))) +
      theme_void(base_size = 11) +
      theme(legend.position = 'right',
            legend.title = element_blank()) +
      coord_polar(clip = 'off', theta = 'x')
    return(p4)
  }
  else
    return(p3)
