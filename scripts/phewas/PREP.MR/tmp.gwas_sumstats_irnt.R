library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)


f.fields = 'data/MR/MR_pheno_ls_UKB.RData'
f.Neale = 'data/MR/UKBB_GWAS_Imputed_v3_201807.tsv'
f.BIG40 = 'data/MR/BIG40_GWAS.csv'
f.phecount = 'data/phe_count.rds'
d.output = 'data/MR/MR_sumstats'


# Load data ---------------------------------------------------------------

# Neale gwas
neale.gwasls = fread(f.Neale)

# BIG40 IDP gwas
big40.gwasls = read.csv(f.BIG40,skip = 1,stringsAsFactors = F) %>% 
  select(-X,-X1) %>% 
  mutate(Pheno=sprintf("%04d", Pheno))

# Pheno-o-i
# fields.single: download
# fields.multiple: download, meta-analyse
# fields.localgwas: no need to download, local gwas
load(f.fields)

# Download single GWAS ----------------------------------------------------
# phenotypes
ls.single.phewas.1 = fields.single$field_used %>%
  paste0('^',.,'$',collapse = '|')
ls.single.phewas.2 = fields.single$field_used %>%
  paste0('^',.,'_',collapse = '|')
ls.single.phewas = c(ls.single.phewas.1,ls.single.phewas.2) %>%
  paste0(.,collapse = '|')
rm(ls.single.phewas.1,ls.single.phewas.2)


# Neale lab
ls.single.neale = neale.gwasls %>%
  .[grepl(ls.single.phewas,.$`Phenotype Code`),] %>%
  filter(Sex=='both_sexes') %>%
  .[grepl('_irnt',.$`Phenotype Code`)] %>% 
  mutate(field_id = gsub('_irnt','',`Phenotype Code`))


# Download GWAS for meta-analysis -----------------------------------------

ls.multiple.phewas.1=c(fields.multiple$field_used_1,fields.multiple$field_used_2) %>% 
  paste0('^',.,'$',collapse = '|')
ls.multiple.phewas.2=c(fields.multiple$field_used_1,fields.multiple$field_used_2) %>% 
  paste0('^',.,'_',collapse = '|')
ls.multiple.phewas = c(ls.multiple.phewas.1,ls.multiple.phewas.2) %>% 
  paste0(.,collapse = '|')
rm(ls.multiple.phewas.1,ls.multiple.phewas.2)

# Neale lab
ls.multiple.neale = neale.gwasls %>% 
  .[grepl(ls.multiple.phewas,.$`Phenotype Code`),] %>% 
  filter(Sex=='both_sexes') %>% 
  .[grepl('_irnt',.$`Phenotype Code`)] %>% 
  mutate(field_id = gsub('_irnt','',`Phenotype Code`))


# Final list of irnt phenotypes -------------------------------------------

ls.irnt = c(ls.single.neale$field_id,ls.multiple.neale$field_id) %>% 
  unique %>% 
  data.frame(field_id=.)

saveRDS(ls.irnt,file='data/MR/ls.irnt.rds')


# Remove non-irnt phenotypes in the MR intermediate folder ----------------

ls.irnt = readRDS('data/MR/ls.irnt.rds')

# remove files
dir.target = 'data/MR/MR_InterFiles/'

ls.f = list.files(path=dir.target,recursive = T,full.names = T) 

irnt.kw = ls.irnt$field_id %>%
  paste0('f\\.',.,'\\.',collapse = '|')

ls.f.to_remove = ls.f %>% 
  .[grepl(irnt.kw,.)]

sapply(ls.f.to_remove,unlink)

# remove info record
instru.info = read_tsv('data/MR/MR_InterFiles/pheno_exposure_outcome/Instrument_info.tsv',col_names=F)

system('mv data/MR/MR_InterFiles/pheno_exposure_outcome/Instrument_info.tsv data/MR/MR_InterFiles/pheno_exposure_outcome/bakup_Instrument_info.tsv')

irnt.kw = ls.irnt$field_id %>%
  paste0('f\\.',.,collapse = '|')
new.instru.info = instru.info %>%   
  .[!grepl(irnt.kw,.$X1),]

write_tsv(new.instru.info,'data/MR/MR_InterFiles/pheno_exposure_outcome/Instrument_info.tsv',col_names=F)



