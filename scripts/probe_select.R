bad_prom <- stratification %>% filter(bad_prom == 'difficult') %>% select(probes)
satel <- stratification %>% filter(satel == 'difficult') %>% select(probes)

low_cg <- cg_density %>% filter(cg == '< 15') %>% select(probes)
high_cg <- cg_density %>% filter(cg == '> 85') %>% select(probes)


tmp <- union_all(bad_prom, satel)
tmp <- union_all(tmp, low_cg)
tmp <- union_all(tmp, high_cg)


df_SWAN <- df_SWAN %>% mutate(annot = ifelse(probes %in% tmp$probes, 'LOW', 'HIGH'))
low <- df_SWAN %>% filter(annot == 'LOW') 
high <- df_SWAN %>% filter(annot == 'HIGH') 


df_SWAN %>% filter(annot == 'LOW') %>% 
  ggplot() +
  geom_point(aes(x = beta_epic, y = beta_np), alpha = 0.2)

cor(low$beta_epic, low$beta_np)
cor(high$beta_epic, high$beta_np)
cor(df_SWAN$beta_epic, df_SWAN$beta_np)




pb <- fread('/Users/lucreziavaleriani/Dropbox/projects/methylation/signal_meth/data/probes/filter_probes.txt', header = FALSE)
pb_df <- df_SWAN %>% filter(probes %in% pb$V1) 

low <- pb_df %>% filter(abs_diff > 0.3 | annot == 'LOW') 
cor(low$beta_epic, low$beta_np)

high <- df_SWAN %>% filter(!(probes %in% low$probes))
cor(high$beta_epic, high$beta_np)




dmr <- fread('/Users/lucreziavaleriani/Desktop/orfeo_LTS/methylation/prova_GEL/LR-Methyation/data/cpg_DM.csv')
dmr <- df_SWAN %>% filter(probes %in% dmr$IlmnID)
dmr %>% filter(abs_diff > 0.3)

cor(dmr$beta_epic, dmr$beta_np)


df_SWAN %>% filter(abs_diff > 0.3)
