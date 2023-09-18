library(optparse)
library(dplyr)
library(tidyr)
library(data.table)
library(patchwork)
library(ggplot2)
library(ggseqlogo)

setwd('./')
source('./scripts/utils.R')

options(bitmapType='cairo')
option_list <- list(
  make_option(c("-l", "--samplelist"), type="character", default=NULL,
              help = "list of samples", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help = "output directory", metavar="character"))

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

out_dir <- file.path(opt$output, opt$sample)
dir.create(out_dir, showWarnings = FALSE)


sample_list <- read.csv('./data/sample_list',
               header = FALSE,
               col.names = c('f_nanopore', 'f_epic', 'sample_name'))

df <- dplyr::tibble()
for (s in seq(1, nrow(sample_list))){
  sample <- sample_list[s,]
  tmp <- create_df(sample$f_nanopore, sample$f_epic)
  id <- rep(sample$sample_name, nrow(tmp))
  tmp <- bind_cols(tmp, sample_ID = id)
  tmp <- tmp %>% dplyr::mutate(abs_diff = beta_epic - beta_np)
  df <- bind_rows(df, tmp) %>% dplyr::select(-score, -canon, -mod, -nread, 
                                             -nread_minus, -canon_minus, -mod_minus, 
                                             -score_minus)
}

intra <- df %>% group_by(probes) %>% 
  mutate(mean_pb = mean(abs_diff)) %>% 
  mutate(median_pb = median(abs_diff)) %>% 
  mutate(median_cov_pb = median(cov))  %>% 
  mutate(nsample = length(unique(sample_ID)))
           
saveRDS(object = intra, file = paste0(out_dir, '/intra_probes.RDS'))




