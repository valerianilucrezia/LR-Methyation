setwd('./')
if (!require("ggseqlogo", quietly = TRUE))
  install.packages("ggseqlogo")

library(optparse)
library(dplyr)
library(data.table)
library(patchwork)
library(ggplot2)
library(ggseqlogo)
options(bitmapType='cairo')

my_palette <- c('palegreen4', 'cadetblue', 'burlywood3', 'darkseagreen', 'darksalmon', 
                'lightpink2', 'indianred','ivory3', 'skyblue4', 'plum4', 
                'tan3', 'lightgoldenrod2', 'firebrick4', 'darkgoldenrod3', 'slateblue4',
                'slategray4', 'darkolivegreen4', 'chocolate4', 'hotpink3', 'steelblue2')


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
kmer <- readRDS('./data/kmer.RDS') 

cs1 <-  ggseqlogo::make_col_scheme(chars=c('A', 'T', 'C', 'G') , 
                                   cols=c('seagreen4', 'brown', 'dodgerblue4', 'goldenrod'))

by <- join_by(probes)
df_kmer <- left_join(df, kmer, by)


get_corr <- function(v1, v2){
  test <- cor.test(v1, v2)
  return(list(test$estimate, test$p.value))
}


analyze_kmer <- function(all, startk, endk, kk){
  mer <- lapply(all$kmer, FUN = function(kmer){
    k <- substr(kmer, 101 - startk, 102 + endk)
  })
  mer <- mer %>% unlist()
  all <- bind_cols(all, n_kmer = mer)
  colnames(all)[ncol(all)] <- kk
  
  c_EN <- c()
  p_EN <- c()
  
  c_PM <- c()
  p_PM <- c()
  
  nprobes <- c()
  
  to_use <- all %>% 
    dplyr::group_by(all[[kk]]) %>% 
    summarize(count = length(unique(probes))) %>% 
    dplyr::filter(count > 2) %>% 
    dplyr::select(`all[[kk]]`) %>% 
    unlist()
  
  for (k in to_use){
    df <- all %>% dplyr::filter(all[[kk]] == k) 
    
    t_EN <- get_corr(df$beta_epic, df$beta_np)
    c_EN <- c(c_EN, t_EN[[1]])
    p_EN <- c(p_EN, t_EN[[2]])
    
    t_PM <- get_corr(df$beta_plus, df$beta_minus)
    c_PM <- c(c_PM, t_PM[[1]])
    p_PM <- c(p_PM, t_PM[[2]])
    
    nprobes <- c(nprobes, nrow(df))
    
  }
  df <- tibble(kmer =  to_use,  
               c_EN = c_EN, 
               p_EN = p_EN, 
               c_PM = c_PM, 
               p_PM = p_PM, 
               nprobes = nprobes)
  
  return(list(all, df))}


res <- analyze_kmer(df_kmer, 2, 2, 'kmer_6')
res_df <- res[[1]]
res_kmer <- res[[2]]
EN_low_kmer <- res_kmer %>% filter(p_EN < 0.01, c_EN < 0.95)
PM_low_kmer <- res_kmer %>% filter(p_PM < 0.01, c_PM < 0.9)

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
  filter(kmer_6 %in% EN_low_kmer$kmer) %>% 
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


plots[['boxplot_PM']] <- res_df %>% filter(kmer_6 %in% PM_low_kmer$kmer) %>% 
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


plots[['barplot_PM']] <- PM_low_kmer %>% ggplot() +
  geom_col(aes(y = kmer, x = nprobes, fill = c_EN)) +
  scale_fill_gradientn(name = 'cor(EPIC, Nanopore)', colors=c("darkslategray","darkslategray4", 'snow3')) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y = element_text(size = 8))

plots[['barplot_EN']]  <- EN_low_kmer %>% ggplot() +
  geom_col(aes(y = kmer, x = nprobes, fill = c_EN)) +
  scale_fill_gradientn(name = 'cor(EPIC, Nanopore)', colors=c("darkslategray","darkslategray4", 'snow3')) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y = element_text(size = 8))

get_bins <- function(df, col1, col2, nbins){
  dx <- 1 / nbins
  df_round <- df %>% 
    dplyr::mutate(col1 = ceiling(df[[`col1`]] / dx), col2 = ceiling(df[[`col2`]] / dx)) %>% 
    dplyr::group_by(col1, col2) %>% 
    dplyr::summarise(n = n())
  return(df_round)
}

nbins <- 120
PM_df <- res_df %>% filter(kmer_6 %in% PM_low_kmer$kmer) 
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


EN_df <- res_df %>% filter(kmer_6 %in% EN_low_kmer$kmer) 
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
  plot_layout(design = layout)

ggsave(paste0(out_dir, '/kmer_', opt$sample, '.pdf'), dpi = 100, height = 14, width = 14, units='in')




