library(optparse)
library(dplyr)
library(data.table)
library(patchwork)
library(ggplot2)
options(bitmapType='cairo')

my_theme <- theme_bw() + theme(
  legend.text = element_text(size=12), 
  title = element_text(size=12, color = 'gray20'),
  axis.title.x = element_text(size = 12, color = 'gray20'),
  axis.text.x = element_text(size = 12, color = 'gray20'),
  axis.title.y = element_text(size = 12, color = 'gray20'),
  axis.text.y = element_text(size = 12, color = 'gray20'))

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
annotation <- annotation %>% 
  dplyr::select(-seqnames, -start, -end, -width, -strand)
by <- join_by(probes)
df_annotation <- left_join(df, annotation, by)

plots <- list()
colors_genes <- list(promoter = 'darkolivegreen',  exon = 'palevioletred3' , 
                     intron = 'orange2', intergenic = 'mediumpurple4')

genes <- df_annotation %>% filter(annot.type == 'genes')
plots[['genes_boxplot']] <- genes %>% 
  ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic - beta_np, fill = annot.id), outlier.shape = NA) +
  geom_hline(aes(yintercept = mean(beta_epic - beta_np)), color = 'ivory4') +
  ylim(-0.5, 0.5) +
  xlab('') +
  ylab('EPIC - Nanopore') + 
  scale_fill_manual(values = colors_genes) +
  my_theme +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


plots[['genes_pm']] <- genes %>% 
  tidyr::pivot_longer(c(beta_np, beta_minus, beta_plus)) %>% 
  ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic -  value, fill = name), outlier.shape = NA) +
  scale_fill_manual(values = c('deepskyblue4', 'darkgoldenrod', 'coral3'), 
                    labels = c('Beta +', 'Beta', 'Beta -'), 
                    name = 'Strand') +
  ylab('EPIC - Nanopore') +
  ylim(-0.5, 0.5) +
  xlab('') +
  my_theme +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

corrs <- c()
nprobes <- c()
for (ids in unique(genes$annot.id)){
  tmp <- genes %>% 
        dplyr::filter(annot.id == ids)
  corrs <- c(corrs, round(cor(tmp$beta_epic, tmp$beta_np),3))
  nprobes <- c(nprobes, length(unique(tmp$probes)))
}
df_genes <- dplyr::tibble(annot = unique(genes$annot.id), 
                  nprobes = nprobes,
                  corr = corrs)




colors_cpgs <- list(shore = 'salmon3', shelf = 'deepskyblue4', 
                    island = 'aquamarine4', inter = 'hotpink4')
cpgs <- df_annotation %>% filter(annot.type == 'cpg')


plots[['cpgs_boxplot']] <- cpgs %>% 
  ggplot() +
  geom_boxplot(aes(x = annot.id, y = beta_epic - beta_np, fill = annot.id), outlier.shape = NA) +
  geom_hline(aes(yintercept = mean(beta_epic - beta_np)), color = 'ivory4') +
  ylim(-0.5, 0.5) +
  xlab('') +
  ylab('EPIC - Nanopore') + 
  scale_fill_manual(name = 'Annotation', values = colors_cpgs) +
  my_theme +
  theme(legend.position = 'none')


corrs <- c()
nprobes <- c()
for (ids in unique(cpgs$annot.id)){
  tmp <- cpgs %>% 
        dplyr::filter(annot.id == ids)
  corrs <- c(corrs, round(cor(tmp$beta_epic, tmp$beta_np),3))
  nprobes <- c(nprobes, length(unique(tmp$probes)))
}
df_cpgs <- dplyr::tibble(annot = unique(cpgs$annot.id), 
                  nprobes = nprobes,
                  corr = corrs)

plots[['cpgs_pm']] <- cpgs %>%  
    tidyr::pivot_longer(c(beta_np, beta_minus, beta_plus)) %>% 
    ggplot() +
    geom_boxplot(aes(x = annot.id, y = beta_epic -  value, fill = name), outlier.shape = NA) +
    geom_hline(aes( yintercept = mean(df$beta_epic - df$beta_np)), color = 'ivory4') + 
    scale_fill_manual(values = c('deepskyblue4', 'darkgoldenrod', 'coral3'), 
                      labels = c('Beta +', 'Beta', 'Beta -'), 
                      name = 'Strand') +
    ylab('EPIC - Nanopore') +
    ylim(-0.5, 0.5) +
    xlab('') +
    my_theme +
    theme(legend.position = 'none')



# Stratification ####
stratification <- readRDS('./data/stratification.RDS')
stratification <- stratification %>% dplyr::select(-chrom, -start, -end)

df_strat <- left_join(df, stratification, by)
plots[['strat_boxplot']] <- df_strat %>% 
  tidyr::pivot_longer(cols = c(chain,seg_dup,bad_prom,mappability,satel,tandem,other,ALL)) %>% 
  ggplot() +
  geom_boxplot(aes(x = name, y= beta_epic - beta_np, fill = value), outlier.shape = NA) +
  geom_hline(aes( yintercept = mean(beta_epic - beta_np)), color = 'ivory4') + 
  scale_fill_manual(values = c('seagreen', 'palevioletred3')) +
  ylim(-0.5, 0.5) +
  ylab('EPIC - Nanopore') +
  xlab('') +
  my_theme +
  theme(legend.position = 'bottom', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


plots[['strat_pm']] <- df_strat %>%  
  tidyr::pivot_longer(c(beta_np, beta_minus, beta_plus)) %>% 
  dplyr::rename(strand = name, beta = value) %>% 
  tidyr::pivot_longer(cols = c(chain,seg_dup,bad_prom,mappability,satel,tandem,other,ALL)) %>% 
  dplyr::filter(value == 'difficult') %>% 
  ggplot() +
  geom_boxplot(aes(x = name, y = beta_epic -  beta, fill = strand), outlier.shape = NA) + 
  geom_hline(aes( yintercept = mean(df$beta_epic - df$beta_np)), color = 'ivory4') + 
  scale_fill_manual(values = c('deepskyblue4', 'darkgoldenrod', 'coral3'), 
                    labels = c('Beta +', 'Beta', 'Beta -'), 
                    name = 'Strand') +
  ylab('EPIC - Nanopore') +
  ggtitle('Difficult') +
  ylim(-0.5, 0.5) +
  xlab('') +
  my_theme +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 



corrs_diff <- c()
corrs_notdiff <- c()
nprobes_diff <- c()
nprobes_notdiff <- c()
for (ids in unique(colnames(df_strat)[21:ncol(df_strat)])){
  tmp_diff <- df_strat %>% 
    dplyr::filter(df_strat[[`ids`]] == 'difficult')
  tmp_not_diff <- df_strat %>% 
    dplyr::filter(df_strat[[`ids`]] == 'not_difficult')
  
  corrs_diff <- c(corrs_diff, round(cor(tmp_diff$beta_epic, tmp_diff$beta_np),3))
  nprobes_diff <- c(nprobes_diff, length(unique(tmp_diff$probes)))
  
  corrs_notdiff <- c(corrs_notdiff, round(cor(tmp_not_diff$beta_epic, tmp_not_diff$beta_np),3))
  nprobes_notdiff <- c(nprobes_notdiff, length(unique(tmp_not_diff$probes)))
}
df_corr_strat <- dplyr::tibble(strat = colnames(df_strat)[21:ncol(df_strat)], 
                         nprobes_diff = nprobes_diff,
                         corr_diff = corrs_diff)
                         # nprobes_not_diff = nprobes_notdiff,
                         # corr_not_diff = corrs_notdiff)



colors_cg <- c('palegreen4', 'cadetblue', 'burlywood3', 'darkseagreen', 'darksalmon', 'lightpink2', 'indianred',
                             'ivory3', 'skyblue4', 'plum4', 'tan3', 'lightgoldenrod2')
                             
cg_density <- readRDS('./data/cg_density.RDS')
cg_density <- cg_density %>% dplyr::select(-chrom, -start, -end)
df_cg <- left_join(df, cg_density)
cgs <- unique(cg_density$cg)[unique(cg_density$cg) != '']

plots[['dens_boxplot']] <- df_cg %>% 
  ggplot() +
  geom_boxplot(aes(x = cg, y= beta_epic - beta_np, fill = cg), outlier.shape = NA) +
  geom_hline(aes( yintercept = mean(beta_epic - beta_np)), color = 'ivory4') + 
  scale_fill_manual(values = colors_cg) +
  ylim(-0.5, 0.5) +
  ylab('EPIC - Nanopore') +
  xlab('') +
  my_theme +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


plots[['dens_pm']] <- df_cg %>%  
  tidyr::pivot_longer(c(beta_np, beta_minus, beta_plus)) %>% 
  ggplot() +
  geom_boxplot(aes(x = cg, y = beta_epic -  value, fill = name), outlier.shape = NA) + 
  geom_hline(aes( yintercept = mean(df$beta_epic - df$beta_np)), color = 'ivory4') + 
  scale_fill_manual(values = c('deepskyblue4', 'darkgoldenrod', 'coral3'), 
                    labels = c('Beta +', 'Beta', 'Beta -'), 
                    name = 'Strand') +
  ylab('EPIC - Nanopore') +
  ylim(-0.5, 0.5) +
  xlab('') +
  my_theme +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 



corrs <- c()
nprobes <- c()
for (ids in sort(unique(df_cg$cg))){
  tmp<- df_cg %>% 
    dplyr::filter(cg == ids)
  
  corrs <- c(corrs, round(cor(tmp$beta_epic, tmp$beta_np),3))
  nprobes <- c(nprobes, length(unique(tmp$probes)))
}
df_corr_cg <- dplyr::tibble(density = sort(unique(df_cg$cg)), 
                            nprobes = nprobes,
                            corr = corrs)

layout <- 'ADGGLL
           BEHHMM
           CF#I#N'

plots$genes_boxplot + 
  plots$genes_pm +
  gridExtra::tableGrob(df_genes) + 
  plots$cpgs_boxplot +
  plots$cpgs_pm +
  gridExtra::tableGrob(df_cpgs) + 
  plots$strat_boxplot +
  plots$strat_pm +
  gridExtra::tableGrob(df_corr_strat) + 
  plots$dens_boxplot +
  plots$dens_pm +
  gridExtra::tableGrob(df_corr_cg) + 
  plot_layout(design = layout)

ggsave(paste0(out_dir, '/annotation_', opt$sample, '.pdf'), dpi = 100, height = 14, width = 18, units='in')
