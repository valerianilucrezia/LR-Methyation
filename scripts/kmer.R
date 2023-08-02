library(ggseqlogo)
cs1 <-  make_col_scheme(chars=c('A', 'T', 'C', 'G') , 
                        cols=c('seagreen4', 'brown', 'dodgerblue4', 'goldenrod'))


analyze_kmer <- function(all, startk, endk, kk){
  mer <- lapply(all$kmer, FUN = function(kmer){
    k <- substr(kmer, 101 - startk, 102 + endk)
  })
  mer <- mer %>% unlist()
  all <- bind_cols(all, n_kmer = mer)
  colnames(all)[ncol(all)] <- kk
  
  c_EN <- c()
  p_EN <- c()
  nprobes <- c()
  
  to_use <- all %>% 
    dplyr::group_by(all[[kk]]) %>% 
    summarize(count = length(unique(probes))) %>% 
    dplyr::filter(count > 2) %>% 
    dplyr::select(`all[[kk]]`) %>% 
    unlist()
  
  for (k in to_use){ #unique(all[[kk]])
    df <- all %>% dplyr::filter(all[[kk]] == k) 
    
    t_EN <- get_corr(df$b_epic.x, df$beta_cov)
    c_EN <- c(c_EN, t_EN[[1]])
    p_EN <- c(p_EN, t_EN[[2]])
    nprobes <- c(nprobes, nrow(df))
    
  }
  df <- tibble(kmer =  to_use,  
               c_EN = c_EN, 
               p_EN = p_EN, 
               nprobes = nprobes)
  
  return(list(all, df))}


l <- analyze_kmer(all, 2, 2, 'kmer_6')
all <- l[[1]]
df <- l[[2]]