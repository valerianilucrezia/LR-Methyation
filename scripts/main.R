if (!require("optparse", quietly = TRUE))
  install.packages("optparse")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("data.table", quietly = TRUE))
  install.packages("data.table")

if (!require("patchwork", quietly = TRUE))
  install.packages("patchwork")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("ggseqlogo", quietly = TRUE))
  install.packages("ggseqlogo")

library(optparse)
library(dplyr)
library(tidyr)
library(data.table)
library(patchwork)
library(ggplot2)
library(ggseqlogo)

setwd('./')
source('./scripts/main_analysis.R')
source('./scripts/kmer.R')
source('./scripts/annotation.R')
source('./scripts/DMR.R')
source('./scripts/utils.R')

options(bitmapType='cairo')
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
dir.create(out_dir, showWarnings = FALSE)

f_nanopore <- opt$nanopore
f_epic <- opt$epic

df <- main_analysis(f_nanopore, f_epic, out_dir, opt$sample)
print('main_analysis done')

annotation(df, out_dir, opt$sample)
print('annotation done')

DMR(df, out_dir, opt$sample)
print('DMR analysis done') 

kmer(df, out_dir, opt$sample)
print('kmer analysis done') 




