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
    k <- stringr::substr(kmer, 101 - startk, 102 + endk)
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
