DMR <- function(df, out_dir, sample_name) {
  df <- df %>% dplyr::filter(cov > 20)
  dmr <- data.table::fread('data/cpg_DM.csv') 
  
  df <- df %>% dplyr::mutate(meth = ifelse(probes %in% dmr$IlmnID, 'DM', 'NDM'))
  df_dmr <- df %>% dplyr::filter(meth == 'DM')
  
  plots <- list()
  plots[['bp_all']] <- df %>% ggplot() +
    geom_boxplot(aes(x = meth, y = beta_epic - beta_np, fill = meth), outlier.shape = NA) +
    xlab('') +
    ylab('EPIC - Nanopore') +
    ylim(-0.5, 0.5) +
    scale_fill_manual(values = c('darkseagreen', 'salmon2')) + 
    ggtitle(paste0('DM probes = ', nrow(df_dmr))) +
    my_theme +
    theme(legend.position = 'none')
  
  
  plots[['scatter']] <- df %>% 
    dplyr::filter(meth == 'DM') %>% 
    ggplot() +
    geom_point(aes(x = beta_epic, y = beta_np), color = 'darkseagreen4', alpha = .3, size = .5) +
    geom_abline(intercept = 0, color = 'ivory4') +
    xlab('EPIC') +
    ylab('Nanopore') +
    my_theme
    
  
  plots[['bp_strand_bias']] <- df %>% 
    tidyr::pivot_longer(c(beta_np, beta_minus, beta_plus)) %>% 
    ggplot() +
    geom_boxplot(aes(x = meth, y = beta_epic -  value, fill = name), outlier.shape = NA) +
    scale_fill_manual(values = c('deepskyblue4', 'darkgoldenrod', 'coral3'), 
                      labels = c('Beta +', 'Beta', 'Beta -'), 
                      name = 'Strand') +
    ylab('EPIC - Nanopore') +
    ylim(-0.5, 0.5) +
    xlab('') +
    my_theme +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  stratification <- readRDS('data/stratification.RDS')
  stratification <- stratification %>% 
    dplyr::select(-chrom, -start, -end)
  by <- join_by(probes)
  df_strat <- left_join(df, stratification, by)
  
  plots[['strat_barplot']] <- df_strat %>% 
    tidyr::pivot_longer(cols = c(chain,seg_dup,bad_prom,mappability,satel,tandem,other,ALL)) %>% 
    group_by(name, value, meth) %>% 
    summarize(count = length(unique(probes))) %>% 
    filter(value == 'difficult') %>% 
    ggplot() +
    geom_col(aes(x = name, y = count, fill = meth), position = 'fill') +  #'fil
    scale_fill_manual(values = c('darkseagreen', 'salmon2')) + 
    ggtitle('Difficult to map') +
    my_theme +
    xlab('') +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  
  
  plots[['strat_boxplot']] <- df_strat %>% 
    tidyr::pivot_longer(cols = c(chain, seg_dup, bad_prom, mappability, satel, tandem, other, ALL)) %>% 
    dplyr::filter(value == 'difficult') %>% 
    ggplot() +
    geom_boxplot(aes(x = name, y= beta_epic - beta_np, fill = meth), outlier.shape = NA) +
    geom_hline(aes( yintercept = mean(beta_epic - beta_np)), color = 'ivory4') + 
    scale_fill_manual(values = c('darkseagreen', 'salmon2')) + 
    ylim(-0.5, 0.5) +
    ylab('EPIC - Nanopore') +
    ggtitle('Difficult to map') +
    xlab('') +
    my_theme +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  
                             
  cg_density <- readRDS('data/cg_density.RDS')
  cg_density <- cg_density %>% dplyr::select(-chrom, -start, -end)
  df_cg <- left_join(df, cg_density)
  cgs <- unique(cg_density$cg)[unique(cg_density$cg) != '']
  
  
  plots[['cg_barplot']] <- df_cg %>%   
    dplyr::group_by(cg, meth) %>% 
    dplyr::summarize(count = length(unique(probes))) %>% 
    ggplot() +
    geom_col(aes(x = cg, y =count, fill = meth), position = 'fill') +
    xlab('') + 
    scale_fill_manual(values = c('darkseagreen', 'salmon2')) + 
    my_theme +   
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  
  plots[['cg_boxplot']] <- df_cg %>% 
    ggplot() +
    geom_boxplot(aes(x = cg, y= beta_epic - beta_np, fill = meth), outlier.shape = NA) +
    geom_hline(aes( yintercept = mean(beta_epic - beta_np)), color = 'ivory4') + 
    scale_fill_manual(values = c('darkseagreen', 'salmon2')) + 
    ylim(-0.5, 0.5) +
    ylab('EPIC - Nanopore') +
    xlab('') +
    my_theme +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  annotation <- readRDS('data/annotation.RDS')
  annotation <- annotation %>% 
    dplyr::select(-seqnames, -start, -end, -width, -strand)
  by <- join_by(probes)
  df_annotation <- left_join(df, annotation, by)
  genes <- df_annotation %>% filter(annot.type == 'genes')
  
  plots[['genes_barplot']] <- genes %>%   
    dplyr::group_by(annot.id, meth) %>% 
    dplyr::summarize(count = length(unique(probes))) %>% 
    ggplot() +
    geom_col(aes(x = annot.id, y =count, fill = meth), position = 'fill') +
    xlab('') + 
    scale_fill_manual(values = c('darkseagreen', 'salmon2')) + 
    my_theme +   
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  
  plots[['genes_boxplot']] <- genes %>% 
    ggplot() +
    geom_boxplot(aes(x = annot.id, y= beta_epic - beta_np, fill = meth), outlier.shape = NA) +
    geom_hline(aes( yintercept = mean(beta_epic - beta_np)), color = 'ivory4') + 
    scale_fill_manual(values = c('darkseagreen', 'salmon2')) + 
    ylim(-0.5, 0.5) +
    ylab('EPIC - Nanopore') +
    xlab('') +
    my_theme +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  cpgs <- df_annotation %>% filter(annot.type == 'cpg')
  plots[['cpg_barplot']] <- cpgs %>%   
    dplyr::group_by(annot.id, meth) %>% 
    dplyr::summarize(count = length(unique(probes))) %>% 
    ggplot() +
    geom_col(aes(x = annot.id, y =count, fill = meth), position = 'fill') +
    xlab('') + 
    scale_fill_manual(values = c('darkseagreen', 'salmon2')) + 
    my_theme +   
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  
  plots[['cpg_boxplot']] <- cpgs %>% 
    ggplot() +
    geom_boxplot(aes(x = annot.id, y= beta_epic - beta_np, fill = meth), outlier.shape = NA) +
    geom_hline(aes( yintercept = mean(beta_epic - beta_np)), color = 'ivory4') + 
    scale_fill_manual(values = c('darkseagreen', 'salmon2')) + 
    ylim(-0.5, 0.5) +
    ylab('EPIC - Nanopore') +
    xlab('') +
    my_theme +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  layout <- 'ABM
             ABM
             CDD
             CDD
             EFF
             EFF
             GH#
             GH#
             IL#
             IL#
             ' 
  
  plots$bp_all +
    plots$bp_strand_bias +
    plots$strat_barplot + 
    plots$strat_boxplot +
    plots$cg_barplot + 
    plots$cg_boxplot +
    plots$genes_barplot +
    plots$genes_boxplot +
    plots$cpg_barplot +
    plots$cpg_boxplot +
    plots$scatter +
    patchwork::plot_layout(design = layout, guides = 'collect') +
    patchwork::plot_annotation(title = sample_name)
  ggsave(paste0(out_dir, '/DMR_', sample_name, '.pdf'),  dpi = 100, height = 15, width = 9, units='in' )
}
