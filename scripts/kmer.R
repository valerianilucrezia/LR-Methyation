kmer <- function(df, out_dir, sample_name){

  my_palette <- c('palegreen4', 'cadetblue', 'burlywood3', 'darkseagreen', 'darksalmon', 
                  'lightpink2', 'indianred','ivory3', 'skyblue4', 'plum4', 
                  'tan3', 'lightgoldenrod2', 'firebrick4', 'darkgoldenrod3', 'slateblue4',
                  'slategray4', 'darkolivegreen4', 'chocolate4', 'hotpink3', 'steelblue2')
  
  df <- df %>% dplyr::filter(cov >= 10)
  kmer <- readRDS('data/kmer.RDS') 
  
  cs1 <-  ggseqlogo::make_col_scheme(chars=c('A', 'T', 'C', 'G') , 
                                     cols=c('seagreen4', 'brown', 'dodgerblue4', 'goldenrod'))
  
  by <- join_by(probes)
  df_kmer <- left_join(df, kmer, by)
  
  res <- analyze_kmer(df_kmer, 2, 2, 'kmer_6')
  res_df <- res[[1]]
  res_kmer <- res[[2]]
  EN_low_kmer <- res_kmer %>% dplyr::filter(p_EN < 0.01, c_EN < 0.95)
  PM_low_kmer <- res_kmer %>% dplyr::filter(p_PM < 0.01, c_PM < 0.9)
  
  EN_low_kmer <- EN_low_kmer %>% arrange(c_EN) %>% slice(1:20)
  PM_low_kmer <- PM_low_kmer %>% arrange(c_PM) %>% slice(1:20)
  
  plots <- list()
  plots[['logo_EN']] <- ggplot() +
    ggseqlogo::geom_logo(EN_low_kmer$kmer %>% unlist(), method = 'prob', col_scheme = cs1) +  
    ylab('') +
    xlab('position') +
    my_theme
  
  plots[['logo_PM']] <- ggplot() +
    ggseqlogo::geom_logo(PM_low_kmer$kmer %>% unlist(), method = 'prob', col_scheme = cs1) +  
    ylab('') +
    xlab('position') +
    my_theme
  
  plots[['boxplot_EN']] <- res_df %>% 
    dplyr::filter(kmer_6 %in% EN_low_kmer$kmer) %>% 
    ggplot() +
    geom_boxplot(aes(x = kmer_6, y = beta_epic - beta_np, fill = kmer_6), outlier.shape = NA) +
    geom_hline(aes(yintercept = mean(beta_epic - beta_np)), color = 'ivory4') +
    ylim(-0.5, 0.5) +
    xlab('') +
    ylab('EPIC - Nanopore') + 
    scale_fill_manual(values = my_palette) +
    my_theme +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  plots[['boxplot_PM']] <- res_df %>% 
    dplyr::filter(kmer_6 %in% PM_low_kmer$kmer) %>% 
    tidyr::pivot_longer(c(beta_np, beta_minus, beta_plus)) %>% 
    ggplot() +
    geom_boxplot(aes(x = kmer_6, y = beta_epic - value, fill = name), outlier.shape = NA) +
    geom_hline(aes(yintercept = mean(res_df$beta_epic - res_df$beta_np)), color = 'gray20') +
    ylim(-0.5, 0.5) +
    xlab('') +
    ylab('EPIC - Nanopore') + 
    scale_fill_manual(values = c('deepskyblue4', 'darkgoldenrod', 'coral3'), 
                      labels = c('Beta +', 'Beta', 'Beta -'), 
                      name = 'Strand') +
    my_theme +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  plots[['barplot_PM']] <- PM_low_kmer %>% 
    ggplot() +
    geom_col(aes(y = kmer, x = nprobes, fill = c_EN)) +
    scale_fill_gradientn(name = 'cor(EPIC, Nanopore)', colors=c("darkslategray","darkslategray4", 'snow3')) +
    my_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          axis.text.y = element_text(size = 8))
  
  plots[['barplot_EN']]  <- EN_low_kmer %>% 
    ggplot() +
    geom_col(aes(y = kmer, x = nprobes, fill = c_EN)) +
    scale_fill_gradientn(name = 'cor(EPIC, Nanopore)', colors=c("darkslategray","darkslategray4", 'snow3')) +
    my_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          axis.text.y = element_text(size = 8))
  
  nbins <- 120
  PM_df <- res_df %>% dplyr::filter(kmer_6 %in% PM_low_kmer$kmer) 
  df_bins_PM <- get_bins(PM_df, 'beta_epic', 'beta_np', nbins)
  plots[['scatter_PM']]  <- PM_df %>% 
    ggplot() +
    geom_bin2d(aes(x = beta_epic, y = beta_np), bins = nbins) +
    scale_fill_gradientn(limits = c(1, max(df_bins_PM$n)),
                         colors = c("lightsteelblue1", "lightsteelblue3", 'lightsteelblue4'),
                         values = scales::rescale(c(1, stats::quantile(df_bins_PM$n, .9), max(df_bins_PM$n)))) +
    xlab('EPIC') +
    ylab('Nanopore') +
    my_theme
  
  
  EN_df <- res_df %>% dplyr::filter(kmer_6 %in% EN_low_kmer$kmer) 
  df_bins_EN <- get_bins(EN_df, 'beta_epic', 'beta_np', nbins)
  plots[['scatter_EN']]  <- EN_df %>% 
    ggplot() +
    geom_bin2d(aes(x = beta_epic, y = beta_np), bins = nbins) +
    scale_fill_gradientn(limits = c(1, max(df_bins_EN$n)),
                         colors = c("lightsteelblue1", "lightsteelblue3", 'lightsteelblue4'),
                         values = scales::rescale(c(1, stats::quantile(df_bins_EN$n, .9), max(df_bins_EN$n)))) +
    xlab('EPIC') +
    ylab('Nanopore') +
    my_theme
  
  layout <- '#BB
             #BB
             ACD
             #CD
             FFF
             FFF
             EGH
             #GH'
  
  plots$logo_EN +
    plots$boxplot_EN +
    plots$barplot_EN +
    plots$scatter_EN +
    plots$logo_PM +
    plots$boxplot_PM +
    plots$barplot_PM +
    plots$scatter_PM + 
    patchwork::plot_layout(design = layout) +
    patchwork::plot_annotation(title = sample_name)
  
  ggsave(paste0(out_dir, '/kmer_', sample_name, '.pdf'), dpi = 100, height = 14, width = 14, units='in')
}




