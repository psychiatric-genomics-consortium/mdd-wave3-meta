library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.fields = args[1] 
f.single = args[2]
d.multiple = args[3]
f.localgwas = args[4]
f.Neale.wget = args[5]
f.Neale.n = args[6]
f.BIG40 = args[7]
f.pipeline.input = args[8]
f.n.ref = args[9]

# f.fields = 'data/MR/MR_pheno_ls_UKB.RData'
# f.single = 'data/MR/single.ls.rds'
# d.multiple = 'data/MR/MR_sumstats'
# f.localgwas = 'data/MR/local.gwas.rds'
# f.Neale.wget = 'data/MR/UKBB_GWAS_Imputed_v3_201807.tsv'
# f.Neale.n = 'data/MR/phenotypes.both_sexes.v2.tsv'
# f.BIG40 = 'data/MR/BIG40_GWAS.csv'
# f.pipeline.input = 'data/MR/pheno_gwas_forMR.rds'
# f.n.ref = 'data/MR/pheno_gwas_n.rds'


# Create processing list --------------------------------------------------
load(f.fields)

ls.single = readRDS(f.single) %>% 
  mutate(loading_type = ifelse(grepl('wget',wget_command),'neale','big40')) %>% 
  select(file_loc,field_tag,loading_type)
ls.multiple = list.files(path=d.multiple,pattern = 'mtag_meta.txt.gz',full.names = T) %>%
  data.frame(file_loc=.) %>% 
  mutate(field_tag=file_loc %>% 
           gsub('_mtag_meta.txt.gz','',.) %>% 
           gsub('data/MR/MR_sumstats/','',.)) %>% 
  mutate(loading_type = 'mtag') %>% 
  select(file_loc,field_tag,loading_type)
ls.localgwas = list.files(path=d.multiple,pattern = 'regenie.gz',full.names = T) %>%
  data.frame(file_loc=.) %>% 
  mutate(field_tag=file_loc %>% 
           gsub('.regenie.gz','',.,fixed = T) %>% 
           gsub('data/MR/MR_sumstats/','',.)) %>% 
  mutate(loading_type = 'local') %>% 
  select(file_loc,field_tag,loading_type)

ls.all = ls.single %>% 
  rbind(.,ls.multiple) %>%
  rbind(.,ls.localgwas) %>% 
  mutate(file_loc=gsub('tsv.bgz','tsv.gz',file_loc,fixed = T)) %>% 
  as.data.frame

saveRDS(ls.all,file=f.pipeline.input)


# Make sample size reference ----------------------------------------------

neale.gwasls = fread(f.Neale.wget) %>% 
  rename(phenotype=`Phenotype Code`) %>% 
  filter(phenotype!='N/A') %>%
  .[!grepl('_irnt',.$phenotype),] %>% 
  mutate(phenotype = gsub('_raw','',phenotype)) %>% 
  left_join(.,fread(f.Neale.n),by='phenotype') %>% 
  filter(!is.na(n_non_missing)) %>% 
  select(phenotype,n=n_non_missing)

# BIG40 IDP gwas
big40.gwasls = read.csv(f.BIG40,skip = 1,stringsAsFactors = F) %>% 
  select(-X,-X1) %>% 
  mutate(Pheno=sprintf("%04d", Pheno)) %>% 
  select(phenotype=UKB.ID,n=N.all.)

gwasls.all = rbind(neale.gwasls,big40.gwasls)

saveRDS(gwasls.all,f.n.ref)