library(optparse)
library(dplyr)
library(data.table)
library(patchwork)
library(ggplot2)
options(bitmapType='cairo')

my_theme <- theme_bw() + theme(
  legend.text = element_text(size=15), 
  title = element_text(size=13, color = 'gray20'),
  axis.title.x = element_text(size = 16, color = 'gray20'),
  axis.text.x = element_text(size = 15, color = 'gray20'),
  axis.title.y = element_text(size = 16, color = 'gray20'),
  axis.text.y = element_text(size = 15, color = 'gray20'))

option_list <- list(
  make_option(c("-o", "--output"), type="character", default=NULL,
              help = "output directory", metavar="character"), 
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help = "sample name", metavar="character"));


opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

out_dir <- file.path(opt$output, opt$sample)
df <- readRDS(paste0(out_dir,'/', opt$sample, '.RDS'))
df <- df %>% filter(cov > 20)

annotation <- readRDS('./data/annotation.RDS')
annotation <- annotation %>% dplyr::select(-seqnames, -start, -end, -width, -strand)
by <- join_by(probes)
df_annotation <- left_join(df, annotation, by)

plots <- list()
colors_genes <- list(promoter = 'darkolivegreen',  exon = 'palevioletred3' , 
                     intron = 'orange2', intergenic = 'mediumpurple4')

genes <- df_annotation %>% filter(annot.type == 'genes')
plots[['genes_boxplot']] <- genes %>% 
  ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic - beta_np, fill = annot.id), outlier.shape = NA) +
  ylim(0.5,-0.5) + 
  xlab('') +
  ylab('EPIC - Nanopore') + 
  scale_fill_manual(values = colors_genes) +
  my_theme +
  theme(legend.position = 'none')


plots[['genes_pm']] <- genes %>% 
  tidyr::pivot_longer(c(beta_np, beta_minus, beta_plus)) %>% 
  ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic -  value, fill = name), outlier.shape = NA) +
  scale_fill_manual(values = c('deepskyblue4', 'darkgoldenrod', 'coral3'), 
                    labels = c('Beta +', 'Beta cov', 'Beta -'), 
                    name = 'Strand') +
  ylab('EPIC - Nanopore') +
  ylim(0.5,-0.5) + 
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
      ggtitle(paste0('Genes_', ids, ': ', round(cor(df$beta_epic, df$beta_np),3), ', ', length(unique(df$probes)))) + 
      xlab('EPIC') +
      ylab('Nanopore') +
      my_theme
})
plots[['genes_scatters']] <- ggpubr::ggarrange(plotlist = genes_scatters, nrow = 2, ncol = 2)


colors_cpgs <- list(shore = 'salmon3', shelf = 'deepskyblue4', 
                    island = 'aquamarine4', inter = 'hotpink4')
cpgs <- df_annotation %>% filter(annot.type == 'cpg')


plots[['cpgs_boxplot']] <- cpgs %>% ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic - beta_np, fill = annot.id), outlier.shape = NA) +
  ylim(0.5,-0.5) + 
  xlab('') +
  ylab('EPIC - Nanopore') + 
  scale_fill_manual(name = 'Annotation', values = colors_cpgs) +
  my_theme +
  theme(legend.position = 'none')

cpgs_scatters <- lapply(unique(cpgs$annot.id), FUN = function(ids){
  df <- cpgs %>% dplyr::filter(annot.id == ids)
  ggplot(df) +
    geom_point(aes(x = beta_epic, y = beta_np), color = colors_cpgs[[ids]], alpha = .3) +
    geom_abline(intercept = 0, color = 'ivory2', linetype='twodash') + 
    geom_smooth(aes(x = beta_epic, y = beta_np), method = "lm", color = 'ivory2', alpha = .2) +
    ggtitle(paste0('CpG_', ids, ': ', round(cor(df$beta_epic, df$beta_np),3), ', ', nrow(df))) + 
    xlab('EPIC') +
    ylab('Nanopore') +
    my_theme
})
plots[['cpgs_scatters']] <- ggpubr::ggarrange(plotlist = cpgs_scatters)

plots[['cpgs_pm']] <- cpgs %>%  tidyr::pivot_longer(c(beta_np, beta_minus, beta_plus)) %>% 
  ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic -  value, fill = name), outlier.shape = NA) +
  scale_fill_manual(values = c('deepskyblue4', 'darkgoldenrod', 'coral3'), 
                    labels = c('Beta +', 'Beta cov', 'Beta -'), 
                    name = 'Strand') +
  ylab('EPIC - Nanopore') +
  ylim(0.5,-0.5) + 
  xlab('') +
  my_theme +
  theme(legend.position = 'bottom')

layout <- 'ABBBB
           #BBBB
           CDDDD
           #DDDD'


plots$genes_boxplot +
  plots$genes_scatters +
  # plots$genes_pm +
  plots$cpgs_boxplot +
  plots$cpgs_scatters +
  # plots$cpgs_pm +
  plot_layout(design = layout)

ggsave(paste0(out_dir, '/annotation_', opt$sample, '.pdf'), dpi = 300, height = 15, width = 12, units='in')
