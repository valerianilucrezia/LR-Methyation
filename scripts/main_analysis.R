main_analysis <- function(f_nanopore, f_epic, out_dir, sample_name){
  np <- data.table::fread(f_nanopore)
  colnames(np) <- c("chrom", "start", "end", "name", "score", "strand", "tstart", "tend", "color", "coverage", "freq", "canon", "mod", "filt")
  
  np <- np %>% dplyr::filter(chrom %in% chrs) %>% 
    dplyr::mutate(beta = mod/(mod+canon)) %>% 
    dplyr::mutate(nread = mod+canon) %>% 
    dplyr::select(-name, -freq, -filt, -tstart, -tend, -color)


  color_strand <- list(plus = 'deepskyblue4',
                       minus = 'coral3',
                       mean = 'darkgoldenrod')
  l_strand <- list(minus = '-', plus = '+')

  plots <- list()
  i <- 1
  pl_strand <- list()
  for (st in c('minus', 'plus')){
    
    df <- np %>% dplyr::filter(strand == l_strand[[st]])  
    pl1 <- df %>% 
      ggplot() +
      geom_histogram(aes(x = nread), fill = color_strand[[st]], alpha = 0.8, bins = 100) +
      xlab('Coverage') + 
      xlim(0, 200) +
      my_theme 
      
    pl2 <- df %>% 
      ggplot() +
      geom_histogram(aes(x = beta), fill = color_strand[[st]], alpha = 0.8, bins = 100) +
      xlab(paste0('Beta', l_strand[[st]])) + 
      xlim(-0.01,1.01) +
      my_theme 
    
    pl_strand[[i]] <- pl1
    i <- i+1
    pl_strand[[i]] <- pl2
    i <- i+1
  }

  np_cov <- np[1:nrow(np)-1] %>% 
    dplyr::mutate(beta_minus = np$beta[2:nrow(np)]) %>% 
    dplyr::mutate(coverage_minus = np$coverage[2:nrow(np)]) %>% 
    dplyr::mutate(nread_minus = np$nread[2:nrow(np)]) %>% 
    dplyr::mutate(canon_minus = np$canon[2:nrow(np)]) %>% 
    dplyr::mutate(mod_minus = np$mod[2:nrow(np)]) %>% 
    dplyr::mutate(score_minus = np$score[2:nrow(np)]) %>% 
    dplyr::rename(beta_plus = beta) %>% 
    dplyr::filter(strand == '+') %>% 
    dplyr::mutate(beta = ((mod+mod_minus)/(mod+mod_minus+canon+canon_minus))) %>% 
    dplyr::mutate(cov =(mod+mod_minus+canon+canon_minus)) %>% 
    dplyr::select(-coverage, -strand, -coverage_minus)
  
  pl1 <- np_cov %>% 
      ggplot() +
      geom_histogram(aes(x = cov), fill = color_strand$mean, alpha = 0.8, bins = 200) +
      xlab('Coverage') + 
      xlim(0, 200) +
      my_theme
    
  pl2 <- np_cov %>% 
      ggplot() +
      geom_histogram(aes(x = beta), fill = color_strand$mean, alpha = 0.8, bins = 200) +
      xlab('Beta') + 
      xlim(-0.01,1.01) +
      my_theme 
  
  pl_strand[[i]] <- pl1
  pl_strand[[i+1]] <- pl2
  pl_strand[[1]] + 
    pl_strand[[2]] +
    pl_strand[[3]] + 
    pl_strand[[4]] +
    pl_strand[[5]] +
    pl_strand[[6]] +
    patchwork::plot_layout(nrow=3, ncol=2)
  ggsave(paste0(out_dir, '/nanopore_', sample_name, '.pdf'), dpi = 300, height=7, width=10, units='in')

  # EPIC analysis ####
  plots <- list()
  
<<<<<<< HEAD
  manifest <- 'data/manifest_epic.RDS'
=======
  manifest <- './data/manifest_epic.RDS'
>>>>>>> refs/remotes/origin/main
  manifest_epic <- readRDS(manifest)
  
  epic <- data.table::fread(f_epic, header = FALSE, skip = 1)
  colnames(epic) <- c('probes', 'beta_epic')
  by_probes <- join_by(probes)
  epic <-  left_join(manifest_epic, epic, by_probes) %>% 
    dplyr::filter(!(is.na(beta_epic)))
  
  # Join and analyze data ####
  epic <- epic %>% 
    dplyr::mutate(id = paste(chrom, start, end, sep = ':')) %>% 
    dplyr::select(-chrom, -start, -end)
  np_cov <- np_cov %>%   # takes time
    dplyr::mutate(id = paste(chrom, start, end, sep = ':')) %>% 
    dplyr::rename(beta_np = beta) 
  by <- join_by(id)
  df <- left_join(epic, np_cov, by)
  df <- df %>% 
    dplyr::filter(!is.na(beta_epic)) %>% 
    dplyr::filter(!(is.na(beta_np))) %>% 
    dplyr::mutate(diff = beta_epic - beta_np,
                  abs_diff = abs(beta_epic - beta_np))
  
  plots[['density']] <- df %>% 
    tidyr::pivot_longer(cols = c(beta_epic, beta_np)) %>% 
    ggplot() +
    geom_density(aes(value, color = name)) +
    scale_color_manual(values = c('turquoise4', 'darkgoldenrod'), 
                       labels = c('EPIC', 'Nanopore'), 
                       name = '') + 
    my_theme + 
    xlim(-0.05, 1.05) + 
    theme(legend.position = 'bottom')
  
  
  plots[['coverage']] <- df %>% 
    ggplot() +
    geom_histogram(aes(x = cov), binwidth = 2, fill = 'darkgoldenrod', alpha = .7) +
    geom_vline(aes(xintercept = mean(df$cov)), color = 'darkgoldenrod1') +
    xlab('Coverage') +
    xlim(0, 150) +
    my_theme
  
  
  corr <- c()
  nprobes <- c()
  diff_mean <- c()
  diff_median <- c()
  for (cv in seq(1, 120)){
    df1 <- df %>% dplyr::filter(cov <= cv)
    df2 <- df %>% dplyr::filter(cov >= cv)
    
    crr <- cor(df1$beta_epic, df1$beta_np)
    corr <- c(corr, crr)
    nprobes <- c(nprobes, length(unique(df2$probes)))
    diff_mean <- c(diff_mean, mean(df1$beta_epic - df1$beta_np))
    diff_median <- c(diff_median, median(df1$beta_epic - df1$beta_np))
  }
  
  cov_cor <-dplyr::tibble(cov = seq(1, 120), 
                    corr = corr, 
                    nprobes = nprobes, 
                    diff_mean = diff_mean,
                    diff_median = diff_median)
  
  plots[['corr1']] <- cov_cor %>% 
    ggplot() +
    geom_point(aes(x = cov, y = corr), color = 'lightsalmon3', size = 0.5) +
    geom_hline(aes(yintercept = cor(df$beta_epic, df$beta_np)), color = 'lightsalmon') +
    geom_vline(aes(xintercept = mean(df$cov)), color = 'lightsalmon4') +
    ylim(0.5,1) +
    ylab('corr(EPIC, Nanopore)') +
    xlab('Coverage') +
    my_theme
  
  
  plots[['corr2']]  <- cov_cor %>% 
    ggplot() +
    geom_point(aes(x = cov, y = nprobes/nrow(df)), color = 'lightsalmon3', size = 0.5) +
    geom_vline(aes(xintercept = mean(df$cov)), color = 'lightsalmon4') +
    ylab('n_probes/tot_probes') +
    xlab('Coverage') +
    xlim(0,120) +
    my_theme
  
  
  df <- df %>% dplyr::filter(cov > 20)
  
  plots[['beta_diff']]  <- df %>% ggplot() +
    geom_histogram(aes(x = beta_epic - beta_np), fill = 'lightsalmon3', alpha = 0.8, binwidth = 0.01) +
    geom_vline(aes(xintercept = mean( beta_epic - beta_np)), color = 'lightsalmon4') +
    xlab('EPIC - Nanopore') +
    my_theme
  
  
  dff <- c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35 ,0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,0.95,  1)
  nprobes <- c()
  for (d in dff){
    tmp <- df %>% dplyr::filter(diff >= d)
    nprobes <- c(nprobes, length(unique(tmp$probes)))
  }
  df_diff <-  dplyr::tibble(diff =  dff,
                    nprobes = nprobes)
  
  plots[['scatter_probes']]  <-  df_diff %>% 
    ggplot() +
    geom_point(aes(x = as.factor(diff), y = nprobes/nrow(df)),color='lightsalmon3' , size = 0.5) +
    xlab('abs(EPIC - Nanopore) <= x') + 
    ylab('n_probes/tot_probes') +
    my_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9)) 
  
  nbins <- 150
  df_bins_en <- get_bins(df, 'beta_epic', 'beta_np', nbins)
  plots[['scatter_en']] <- df %>%
    ggplot() +
    geom_bin2d(aes(x = beta_epic, y = beta_np), bins = nbins) +
    scale_fill_gradientn(limits=c(1, max(df_bins_en$n)), 
                         colors=c("lightyellow2","darkgoldenrod1", 'darkgoldenrod'),
                         values = scales::rescale(c(1, stats::quantile(df_bins_en$n, .9), max(df_bins_en$n)))) +
    geom_abline(intercept = 0, color = 'ivory4', linetype='twodash') +
    geom_smooth(aes(x = beta_epic, y = beta_np), method = "lm", color = 'ivory4', alpha = .2) +
    xlab('EPIC') + 
    ylab('Nanopore') + 
    ggtitle(paste0('Nprobes = ', nrow(df), ' Corr = ', round(cor(df$beta_epic, df$beta_np), 3)) ) + 
    my_theme +
    theme(legend.position = 'none')
  
  
  df_bins_plus <- get_bins(df, 'beta_epic', 'beta_plus', nbins)
  plots[['scatter_plus']] <- df %>%
    ggplot() +
    geom_bin2d(aes(x = beta_epic, y = beta_plus), bins = nbins) +
    scale_fill_gradientn(limits=c(1, max(df_bins_plus$n)), 
                         colors=c("powderblue","deepskyblue3", 'deepskyblue4'),
                         values = scales::rescale(c(1, stats::quantile(df_bins_plus$n, .9), max(df_bins_plus$n)))) +
    geom_abline(intercept = 0, color = 'ivory4', linetype='twodash') +
    geom_smooth(aes(x = beta_epic, y = beta_plus), method = "lm", color = 'ivory4', alpha = .2) +
    xlab('EPIC') + 
    ylab('Beta +') + 
    my_theme +
    theme(legend.position = 'none')
  
  df_bins_minus <- get_bins(df, 'beta_epic', 'beta_minus', nbins)
  plots[['scatter_minus']] <- df %>% 
    ggplot() +
    geom_bin2d(aes(x = beta_epic, y = beta_minus), bins = nbins) +
    scale_fill_gradientn(limits=c(1, max(df_bins_minus$n)), 
                         colors=c("lavenderblush2","coral1", 'coral3'),
                         values =scales::rescale(c(1, stats::quantile(df_bins_minus$n, .9), max(df_bins_minus$n)))) +
    geom_abline(intercept = 0, color = 'ivory4', linetype='twodash') +
    geom_smooth(aes(x = beta_epic, y = beta_minus), method = "lm", color = 'ivory4', alpha = .2) +
    xlab('EPIC') + 
    ylab('Beta -') + 
    my_theme +
    theme(legend.position = 'none')
  
  
  plots[['coverage']] <- df %>% 
    ggplot() +
    geom_histogram(aes(x = cov), binwidth = 0.5, fill = 'darkgoldenrod', alpha = .7) +
    xlab('Coverage') +
    xlim(0, 150) +
    my_theme
  
  
  df_bins_bias <- get_bins(df, 'beta_plus', 'beta_minus', nbins)
  plots[['scatter_bias']] <- df %>% 
    ggplot() +
    geom_bin2d(aes(x = beta_plus, y = beta_minus), bins = nbins) +
    scale_fill_gradientn(limits = c(1, max(df_bins_bias$n)),
                         colors = c("thistle3", "mediumpurple3", 'mediumpurple4'),
                         values =scales::rescale(c(1, stats::quantile(df_bins_bias$n, .9), max(df_bins_bias$n)))) +
    xlab('Beta+') + 
    ylab('Beta-') + 
    my_theme +
    theme(legend.position = 'none')
  
  plots[['dist_bias']] <- df %>% 
    ggplot() +
    geom_histogram(aes(x = beta_plus - beta_minus), binwidth = 0.02, fill = 'mediumpurple4', alpha = 0.7) +
    geom_vline(aes(xintercept = mean(df$beta_plus - df$beta_minus, na.rm  = TRUE)), color = 'mediumpurple1') +
    xlab('Beta+ - Beta-') +
    xlim(-1.01, 1.01) +
    my_theme +
    theme(legend.position = 'none')
  
  
  df_filter <- df %>% 
    dplyr::mutate(diff_strand = abs(beta_plus - beta_minus)) %>% 
    dplyr::filter(diff_strand > 0.3)
  
  nbins <- 120
  df_bins_bias_beta <- get_bins(df_filter, 'beta_epic', 'beta_np', nbins)
  plots[['scatter_bias_beta']] <- df_filter %>% 
    ggplot() +
    geom_bin2d(aes(x = beta_epic, y = beta_np), bins = nbins) +
    geom_abline(intercept = 0, color = 'ivory4', linetype='twodash') +
    scale_fill_gradientn(limits = c(1, max(df_bins_bias_beta$n)),
                         colors=c("lightyellow2","darkgoldenrod1", 'darkgoldenrod'),
                         values =scales::rescale(c(1, stats::quantile(df_bins_bias_beta$n, .9), max(df_bins_bias_beta$n)))) +
    xlab('EPIC') + 
    ylab('Nanopore') + 
    my_theme +
    theme(legend.position = 'none')
  
  
  
  layout <- "AFE
             BDC
             MGH
             NIL
            "
  
  plots$density +
    plots$coverage +
    plots$corr2 +
    plots$corr1 +
    plots$beta_diff +
    plots$scatter_en +
    plots$scatter_plus +
    plots$scatter_minus +
    plots$scatter_bias +
    plots$dist_bias +
    plots$scatter_probes + 
    plots$scatter_bias_beta + 
    patchwork::plot_layout(design = layout) +
    patchwork::plot_annotation(title = opt$sample)
  
  
  ggsave(paste0(out_dir, '/analysis_', sample_name, '.pdf'), dpi = 200, height = 12, width = 10, units='in')
  return(df)
}



