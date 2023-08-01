option_list <- list(
  make_option(c("-out", "--output"), type="character", default=NULL,
              help = "output directory", metavar="character"), 
  make_option(c("-sample", "--sample"), type="character", default=NULL,
              help = "sample name", metavar="character"));


opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

out_dir <- paste0(opt$out_dir, opt$sample, '/')
df <- readRDS(paste0(opt$out_dir, opt$sample, '.RDS'))

annotation <- readRDS('/orfeo/LTS/LADE/LT_storage/lvaleriani/methylation/prova_GEL/data/annotation.RDS')
annotation <- annotation %>% select(-seqnames, -start, -end, -width, -strand)
by <- join_by(probes)
df_annotation <- left_join(df, annotation, by)


colors_genes <- list(promoter = 'darkolivegreen',  exon = 'palevioletred3' , 
                     intron = 'orange2', intergenic = 'mediumpurple4')

genes <- df_annotation %>% filter(annot.type == 'genes')
genes_boxplot <- genes %>% 
  ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic - beta_np, fill = annot.id), outlier.shape = NA) +
  ylim(0.5,-0.5) + 
  xlab('') +
  ylab('EPIC - Nanopore') + 
  scale_fill_manual(values = colors_genes) +
  my_theme +
  theme(legend.position = 'none')


genes_pm <- genes %>% 
  tidyr::pivot_longer(c(beta_np, beta_minus, beta_plus)) %>% 
  ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic -  value, fill = name)) +
  scale_fill_manual(values = c('deepskyblue4', 'darkgoldenrod', 'coral3'), 
                    labels = c('Beta +', 'Beta cov', 'Beta -'), 
                    name = 'Strand') +
  ylab('EPIC - Nanopore') +
  xlab('') +
  my_theme +
  theme(legend.position = 'bottom')


genes_scatters <- lapply(unique(genes$annot.id), FUN = function(ids){
  df <- genes %>% 
    dplyr::filter(annot.id == ids)
    ggplot(df) +
      geom_point(aes(x = beta_epic, y = beta_np), color = colors_genes[[ids]], alpha = .3) +
      geom_abline(intercept = 0, color = 'ivory2', linetype='twodash') + 
      geom_smooth(aes(x = beta_epic, y = beta_np), method = "lm", color = 'ivory2', alpha = .2) +
      ggtitle(paste0('Genes_', ids, ': correlation = ', round(cor(df$beta_epic, df$beta_np),3), ', nprobes = ', length(unique(df$probes)))) + 
      xlab('EPIC') +
      ylab('Nanopore') +
      my_theme
})
ggpubr::ggarrange(plotlist = genes_scatters, nrow = 2, ncol = 2)


colors_cpgs <- list(shore = 'salmon3', shelf = 'deepskyblue4', 
                    island = 'aquamarine4', inter = 'hotpink4')
cpgs <- df_annotation %>% filter(annot.type == 'cpg')


cpgs_boxplot <- cpgs %>% ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic - beta_np, fill = annot.id), outlier.shape = NA) +
  ylim(0.5,-0.5) + 
  xlab('') +
  ylab('EPIC - NanoporeMean') + 
  scale_fill_manual(name = 'Annotation', values = colors_cpgs) +
  my_theme +
  theme(legend.position = 'none')

cpgs_scatters <- lapply(unique(cpgs$annot.id), FUN = function(ids){
  df <- cpgs %>% dplyr::filter(annot.id == ids)
  ggplot(df) +
    geom_point(aes(x = beta_epic, y = beta_np), color = colors_cpgs[[ids]], alpha = .3) +
    geom_abline(intercept = 0, color = 'ivory2', linetype='twodash') + 
    geom_smooth(aes(x = beta_epic, y = beta_np), method = "lm", color = 'ivory2', alpha = .2) +
    ggtitle(paste0('CpG_', ids, ': correlation = ', round(cor(df$beta_epic, df$beta_np),3), ', nprobes = ', nrow(df))) + 
    xlab('EPIC') +
    ylab('Nanopore') +
    my_theme
})
ggpubr::ggarrange(plotlist = cpgs_scatters)

cpgs_pm <- cpgs %>%  tidyr::pivot_longer(c(beta_np, betas_minus, beta_plus)) %>% 
  ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic -  value, fill = name)) +
  scale_fill_manual(values = c('deepskyblue4', 'darkgoldenrod', 'coral3'), labels = c('Beta +', 'Beta cov', 'Beta -'), name = 'Strand') +
  ylab('EPIC - Nanopore') +
  xlab('') +
  my_theme +
  theme(legend.position = 'bottom')

