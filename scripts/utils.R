my_theme <- theme_bw() + theme(
  legend.text = element_text(size=12), 
  title = element_text(size=12, color = 'gray20'),
  axis.title.x = element_text(size = 12, color = 'gray20'),
  axis.text.x = element_text(size = 12, color = 'gray20'),
  axis.title.y = element_text(size = 12, color = 'gray20'),
  axis.text.y = element_text(size = 12, color = 'gray20'))

chrs <- paste0('chr', 1:22)

get_bins <- function(df, col1, col2, nbins){
  dx <- 1 / nbins
  df_round <- df %>% 
    dplyr::mutate(col1 = ceiling(df[[`col1`]] / dx), col2 = ceiling(df[[`col2`]] / dx)) %>% 
    dplyr::group_by(col1, col2) %>% 
    dplyr::summarise(n = n())
  return(df_round)
}


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
    dplyr::summarize(count = length(unique(probes))) %>% 
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
  df <- dplyr::tibble(kmer =  to_use,  
                      c_EN = c_EN, 
                      p_EN = p_EN, 
                      c_PM = c_PM, 
                      p_PM = p_PM, 
                      nprobes = nprobes)
  
  return(list(all, df))}


create_df <- function(f_nanopore, f_epic){
  np <- data.table::fread(f_nanopore)
  colnames(np) <- c("chrom", "start", "end", "name", "score", "strand", "tstart", "tend", "color", "coverage", "freq", "canon", "mod", "filt")
  
  np <- np %>% dplyr::filter(chrom %in% chrs) %>% 
    dplyr::mutate(beta = mod/(mod+canon)) %>% 
    dplyr::mutate(nread = mod+canon) %>% 
    dplyr::select(-name, -freq, -filt, -tstart, -tend, -color)
  
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
  
  manifest <- 'data/manifest_epic.RDS'
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
  df <- df %>% dplyr::filter(cov >= 10)
  return(df)
}
