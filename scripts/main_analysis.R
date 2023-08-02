library(optparse)
library(dplyr)
library(data.table)
library(patchwork)
library(ggplot2)

setwd('./')
options(bitmapType='cairo')

# Define #### 
my_theme <- theme_bw() + theme(
  legend.text = element_text(size=15), 
  title = element_text(size=13, color = 'gray20'),
  axis.title.x = element_text(size = 16, color = 'gray20'),
  axis.text.x = element_text(size = 15, color = 'gray20'),
  axis.title.y = element_text(size = 16, color = 'gray20'),
  axis.text.y = element_text(size = 15, color = 'gray20'))
chrs <- paste0('chr', 1:22)

option_list <- list(
  make_option(c("-n", "--nanopore"), type="character", default=NULL,
              help = "nanopore file", metavar="character"),
  make_option(c("-e", "--epic"), type="character", default=NULL,
              help = "epic file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help = "output directory", metavar="character"), 
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help = "sample name", metavar="character"));

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

out_dir <- file.path(opt$output, opt$sample)
dir.create(out_dir)

f_nanopore <- opt$nanopore
f_epic <- opt$epic


# Nanopore analysis ####
np <- fread(f_nanopore)
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
  
  df <- np %>% filter(strand == l_strand[[st]])  
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
  mutate(beta_minus = np$beta[2:nrow(np)]) %>% 
  mutate(coverage_minus = np$coverage[2:nrow(np)]) %>% 
  mutate(nread_minus = np$nread[2:nrow(np)]) %>% 
  mutate(canon_minus = np$canon[2:nrow(np)]) %>% 
  mutate(mod_minus = np$mod[2:nrow(np)]) %>% 
  mutate(score_minus = np$score[2:nrow(np)]) %>% 
  rename(beta_plus = beta) %>% 
  filter(strand == '+') %>% 
  mutate(beta = ((mod+mod_minus)/(mod+mod_minus+canon+canon_minus))) %>% 
  mutate(cov =(mod+mod_minus+canon+canon_minus)) %>% 
  select(-coverage, -strand, -coverage_minus)

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
pl_strand[[1]] + pl_strand[[2]] +
  pl_strand[[3]] + pl_strand[[4]] +
  pl_strand[[5]] + pl_strand[[6]] +
  plot_layout(nrow=3, ncol=2)
ggsave(paste0(out_dir, '/nanopore_', opt$sample, '.pdf'), dpi = 300, height=7, width=10, units='in')

# EPIC analysis ####
plots <- list()

manifest <- '../data/manifest_epic.RDS'
manifest_epic <- readRDS(manifest)

epic <- fread(f_epic, header = FALSE, skip = 1)
colnames(epic) <- c('probes', 'beta_epic')
by_probes <- join_by(probes)
epic <-  left_join(manifest_epic, epic, by_probes) %>% 
  filter(!(is.na(beta_epic)))

# Join and analyze data ####
epic <- epic %>% 
  mutate(id = paste(chrom, start, end, sep = ':')) %>% select(-chrom, -start, -end)
np_cov <- np_cov %>%   # takes time
  mutate(id = paste(chrom, start, end, sep = ':')) %>% 
  dplyr::rename(beta_np = beta) 
by <- join_by(id)
df <- left_join(epic, np_cov, by)
df <- df %>% filter(!is.na(beta_epic)) %>% filter(!(is.na(beta_np)))
saveRDS(df, paste0(out_dir,'/', opt$sample, '.RDS'))

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
  xlab('Coverage') +
  xlim(0, 150) +
  my_theme


corr <- c()
nprobes <- c()
diff_mean <- c()
diff_median <- c()
for (cv in seq(1, 120)){
  df1 <- df %>% filter(cov <= cv)
  df2 <- df %>% filter(cov >= cv)
  
  crr <- cor(df1$beta_epic, df1$beta_np)
  corr <- c(corr, crr)
  nprobes <- c(nprobes, length(unique(df2$probes)))
  diff_mean <- c(diff_mean, mean(df1$beta_epic - df1$beta_np))
  diff_median <- c(diff_median, median(df1$beta_epic - df1$beta_np))
}

cov_cor <- tibble(cov = seq(1, 120), 
                  corr = corr, 
                  nprobes = nprobes, 
                  diff_mean = diff_mean,
                  diff_median = diff_median)

plots[['corr1']] <- cov_cor %>% 
  ggplot() +
  geom_point(aes(x = cov, y = corr), color = 'lightsalmon3') +
  geom_hline(aes(yintercept = cor(df$beta_epic, df$beta_np)), color = 'lightsalmon') +
  geom_vline(aes(xintercept = mean(df$cov)), color = 'lightsalmon4') +
  ylim(0.5,1) +
  ylab('corr(EPIC, Nanopore)') +
  xlab('Coverage') +
  my_theme


plots[['corr2']]  <- cov_cor %>% 
  ggplot() +
  geom_point(aes(x = cov, y = nprobes), color = 'lightsalmon3') +
  geom_vline(aes(xintercept = mean(df$cov)), color = 'lightsalmon4') +
  ylab('Nprobes') +
  xlab('Coverage') +
  xlim(0,120) +
  my_theme


df <- df %>% filter(cov > 20)

plots[['beta_diff']]  <- df %>% ggplot() +
  geom_histogram(aes(x = beta_epic - beta_np), fill = 'lightsalmon3', alpha = 0.8, binwidth = 0.01) +
  geom_vline(aes(xintercept = mean( beta_epic - beta_np)), color = 'lightsalmon4') +
  xlab('EPIC - Nanopore') +
  my_theme


plots[['scatter_en']] <- df %>%
  ggplot() +
  geom_bin2d(aes(x = beta_epic, y = beta_np), bins = 150) +
  scale_fill_gradientn(limits=c(1, 2500), colors=c("lightyellow2","darkgoldenrod1", 'darkgoldenrod'),
                       values = scales::rescale(c(1,1000,25000))) +
  geom_abline(intercept = 0, color = 'ivory4', linetype='twodash') +
  geom_smooth(aes(x = beta_epic, y = beta_np), method = "lm", color = 'ivory4', alpha = .2) +
  
  xlab('EPIC') + 
  ylab('Nanopore') + 
  ggtitle(round(cor(df$beta_epic, df$beta_np),3)) +
  my_theme +
  theme(legend.position = 'none')


plots[['scatter_plus']] <- df %>%
  ggplot() +
  geom_bin2d(aes(x = beta_epic, y = beta_plus), bins = 150) +
  scale_fill_gradientn(limits=c(1, 2500), colors=c("powderblue","deepskyblue3", 'deepskyblue4'),
                       values = scales::rescale(c(1,1000,25000))) +
  geom_abline(intercept = 0, color = 'ivory4', linetype='twodash') +
  geom_smooth(aes(x = beta_epic, y = beta_plus), method = "lm", color = 'ivory4', alpha = .2) +
  xlab('EPIC') + 
  ylab('Beta +') + 
  my_theme +
  theme(legend.position = 'none')


plots[['scatter_minus']] <- df %>% 
  ggplot() +
  geom_bin2d(aes(x = beta_epic, y = beta_minus), bins = 150) +
  scale_fill_gradientn(limits=c(1, 2500), 
                       colors=c("lavenderblush2","coral1", 'coral3'),
                       values = scales::rescale(c(1,1000,25000))) +
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

plots[['scatter_bias']] <- df %>% 
  ggplot() +
  geom_bin2d(aes(x = beta_plus, y = beta_minus), bins = 150) +
  scale_fill_gradientn(limits = c(1, 2500),
                       colors = c("thistle3", "mediumpurple3", 'mediumpurple4'),
                       values = scales::rescale(c(1,100,2500))) +
  xlab('Beta+') + 
  ylab('Beta-') + 
  my_theme +
  theme(legend.position = 'none')

plots[['dist_bias']] <- df %>% 
  ggplot() +
  geom_histogram(aes(x = beta_plus - beta_minus), binwidth = 0.05, fill = 'mediumpurple4', alpha = 0.7) +
  xlab('Beta+ - Beta-') +
  xlim(-1.01, 1.01) +
  my_theme +
  theme(legend.position = 'none')

layout <- "ABF
           DCE
           GH#
           IL#
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
  plot_layout(design = layout) 

ggsave(paste0(out_dir, '/analysis_', opt$sample, '.pdf'), dpi = 300, height = 12, width=10, units='in')



