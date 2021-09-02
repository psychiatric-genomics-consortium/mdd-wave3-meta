library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dat_dir = args[1] 
f.dic_dir = args[2]
f.output_data = args[3]

# f.dat_dir = '/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Collab/mdd-meta/results/phewas/phewas_out_Body_MRI.rds'
# f.dic_dir = '/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Collab/mdd-meta/results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt'
# f.category = '/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Collab/mdd-meta/data/phewas_categories.tsv'



# Load data ---------------------------------------------------------------

# Load fields
d.path = f.dic_dir %>% strsplit(.,'/') %>% unlist %>% 
  .[1:(length(.)-1)] %>% 
  paste0(.,collapse = '/')

fields.all = list.files(d.path,pattern='^fields.final',full.names = T,include.dirs = T) %>% 
  as.list %>%
  lapply(fread,stringsAsFactors=F,header=T,sep='\t',quote="") %>%
  lapply(as.data.frame) %>% 
  lapply(mutate,field_used=as.character(field_used),FieldID=as.character(FieldID)) %>% 
  Reduce(function(dtf1,dtf2) bind_rows(dtf1,dtf2[,colnames(dtf1)]), .)

# Load all results
g.path = f.dat_dir %>% strsplit(.,'/') %>% unlist %>% 
  .[1:(length(.)-1)] %>% 
  paste0(.,collapse = '/')

res.targetpT = list.files(g.path,pattern='^phewas_out',full.names = T,include.dirs = T) %>% 
  lapply(readRDS) %>%
  Reduce(function(dtf1,dtf2) bind_rows(dtf1,dtf2), .) %>% 
  mutate(p.corrected = p.adjust(p.value,method='fdr')) %>% 
  mutate(mod_name=as.character(mod_name)) 

# Find target pt: most predictive of MDD_CIDI
tmp.mdd_cidi = res.targetpT %>% 
  .[grepl('^MDD3',.$factor),] %>% 
  .[grep('MDD_CIDI',.$dependent),] %>% 
  .[order(abs(.$beta),decreasing = T),]

target.pT = tmp.mdd_cidi[1,'factor']

res.targetpT = res.targetpT %>% 
  .[grepl(target.pT,.$factor),] %>% 
  .[grepl('^MDD3',.$factor),]

# Annotate results
# ref.category = fread(f.category,header=F,sep='\n') %>% 
#   as.data.frame %>% .$V1

res.targetpT.annot = res.targetpT %>% 
  left_join(.,fields.all,by=c('dependent'='field_tag')) %>% 
  select(category,Field,FieldID,field_used,field_tag=dependent,factor,beta,std,t.value,p.value,p.corrected) %>% 
  filter(p.corrected<0.05) %>% 
  filter(category!='Mental health')

# Parse fields
# fields.single: single fields
# fields.multiple: average of two fields (left + right)
# fields.localgwas: neither of the above, need to do gwas locally (some special items in cognitive tasks)

fields.single = res.targetpT.annot %>% 
  .[.$FieldID==.$field_used,] %>% 
  select(Field,FieldID,field_used)

fields.multiple = res.targetpT.annot %>% 
  .[!.$FieldID %in% fields.single$FieldID,] %>% 
  .[!grepl('.diffs$',.$field_tag),]

ls.multiple = fields.multiple$field_used %>% 
  as.list %>% 
  lapply(strsplit,',') %>% 
  bind_cols %>% 
  t
fields.multiple = fields.multiple %>% 
  mutate(field_used_1=ls.multiple[,1],field_used_2=ls.multiple[,2])

fields.localgwas = res.targetpT.annot %>% 
  .[!.$FieldID %in% fields.single$FieldID,] %>% 
  .[!.$FieldID %in% fields.multiple$FieldID,]


# Save poi list -----------------------------------------------------------
save(fields.single,fields.multiple,fields.localgwas,fields.all,file=f.output_data)