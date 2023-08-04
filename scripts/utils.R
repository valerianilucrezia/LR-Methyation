my_theme <- theme_bw() + theme(
  legend.text = element_text(size=12), 
  title = element_text(size=12, color = 'gray20'),
  axis.title.x = element_text(size = 12, color = 'gray20'),
  axis.text.x = element_text(size = 12, color = 'gray20'),
  axis.title.y = element_text(size = 12, color = 'gray20'),
  axis.text.y = element_text(size = 12, color = 'gray20'))

chrs <- paste0('chr', 1:22)