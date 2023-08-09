## Helper functions
get_ticks <- function(data, col_name) {
  n1 <- floor(log10(range(data[col_name])))
  if (n1[1] == -Inf) n1[1] = 0
  pow <- seq(n1[1], n1[2]+1)
  ticks <- as.vector(sapply(pow, function(p) (c(1,5)*10^p)))
  return(ticks)
}

## plot jaccard rank plot.
plotJaccard = function(df, outD, jaccard_threshold)
{
    target.file = paste0(outD, 'beadsnum_freq.pdf')

    jaccard = sample(df$jaccard, 50000)
    jaccard = jaccard[order(jaccard, decreasing = T)]

    ## jaccard threshold is computed by knee call.

    jaccard.mut = data.frame(jaccard_frag = jaccard) %>%
        arrange(desc(jaccard_frag)) %>% 
        mutate(rank = 1:n(), Merged = jaccard_frag > jaccard_threshold)
    #ticks_at_y = get_ticks(jaccard.mut, 'jaccard_frag')
    ticks_at_x = get_ticks(jaccard.mut, 'rank')

    jaccard_plot = jaccard.mut %>% 
      ggplot(aes(x = rank, y = jaccard_frag)) + 
      #scale_y_log10(breaks = ticks_at_y, labels = as.numeric(ticks_at_y)) +
      scale_x_log10(breaks = ticks_at_x, labels = as.integer(ticks_at_x)) +
      xlab("rank-sorted barcode pairs") +
      ylab("jaccard index") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle=90),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      geom_point(aes(color = Merged), size = 0.6) +
      scale_color_manual(values = c('black', "#3B99D4")) +
      theme(legend.position = c(0.8, 0.8)) +
      labs(color='Merge')
    
    ggsave(jaccard_plot, filename = target.file, width = 6, height = 6, useDingbats = FALSE)
}
