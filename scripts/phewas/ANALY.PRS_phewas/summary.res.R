library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(ggplot2)
library(optparse)


# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dat_dir = args[1]                              
f.dic_dir = args[2]
target.pT = args[3]
f.category = args[4]

#f.output_res = args[4]

# f.dat_dir = 'results/phewas/phewas_out_Body_MRI.rds'
# f.dic_dir = 'results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt'
# target.pT = '0.1'
# f.category = 'data/phewas_categories.tsv'


# Load res and data dictionary --------------------------------------------

target.pT = target.pT %>% 
  gsub('\\.','\\\\.',.) %>% 
  paste0('Pt_',.)

d.path = f.dic_dir %>% strsplit(.,'/') %>% unlist %>% 
  .[1:(length(.)-1)] %>% 
  paste0(.,collapse = '/')

fields.all = list.files(d.path,pattern='^fields.final',full.names = T,include.dirs = T) %>% 
  lapply(fread,stringsAsFactors=F,header=T,sep='\t',quote="\"") %>%
  lapply(as.data.frame) %>% 
  lapply(mutate,field_used=as.character(field_used),FieldID=as.character(FieldID)) %>% 
  Reduce(function(dtf1,dtf2) bind_rows(dtf1,dtf2[,colnames(dtf1)]), .)


g.path = f.dat_dir %>% strsplit(.,'/') %>% unlist %>% 
  .[1:(length(.)-1)] %>% 
  paste0(.,collapse = '/')

res.targetpT = list.files(g.path,pattern='^phewas_out',full.names = T,include.dirs = T) %>% 
  lapply(readRDS) %>%
  Reduce(function(dtf1,dtf2) bind_rows(dtf1,dtf2), .) %>% 
  mutate(p.corrected = p.adjust(p.value,method='fdr')) %>% 
  mutate(mod_name=as.character(mod_name)) %>% 
  .[grepl(target.pT,.$factor),] 

ref.category = fread(f.category,header=F,sep='\n') %>% 
  as.data.frame %>% .$V1

res.targetpT.annot = res.targetpT %>% 
  left_join(.,fields.all,by=c('dependent'='field_tag')) %>% 
  arrange(match(.$category,ref.category)) %>% 
  select(category,Field,field=dependent,factor,beta,std,t.value,p.value,p.corrected)



