library(dplyr)
library(readr)

ls.single.neale = readRDS('data/MR/single.ls.rds') %>% 
  .[!grepl('IDP',.$file_loc)]

ls.multiple.neale = read_tsv('data/MR/mtag_list',col_names=F) %>% 
  select(f1=X1,f2=X2,f.meta=X3) %>% 
  as.data.frame %>% 
  .[!grepl('IDP',.$f1),]

field.single = ls.single.neale %>% .$field_tag

field.multiple = ls.multiple.neale$f.meta %>% 
  strsplit(.,split = '/') %>% 
  lapply(tail,n=1) %>% unlist

field.neale = c(field.single,field.multiple)

saveRDS(field.neale,'data/fields.neale.rds')