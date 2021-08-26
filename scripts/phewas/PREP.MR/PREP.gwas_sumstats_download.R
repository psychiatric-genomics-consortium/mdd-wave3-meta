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
f.Neale = args[2]
f.BIG40 = args[3]
d.output = args[4]

# Neale list: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=178908679
# BIG40 list: https://open.win.ox.ac.uk/ukbiobank/big40/BIG40-IDPs_v4/IDPs.html
# Convert BIG40 to csv: https://www.convertcsv.com/html-table-to-csv.htm

# f.fields = 'data/MR/MR_pheno_ls_UKB.RData'
# f.Neale = 'data/MR/UKBB_GWAS_Imputed_v3_201807.tsv'
# f.BIG40 = 'data/MR/BIG40_GWAS.csv'
# d.output = 'data/MR/MR_sumstats'

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


# Download variant info ---------------------------------------------------

# variants info for Neale's
if (!file.exists('data/MR/MR_sumstats/Neale/variants.tsv.bgz')&
    !file.exists('data/MR/MR_sumstats/Neale/variants.txt')){
  system('wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz')
  system('mv variants.tsv.bgz data/MR/MR_sumstats/Neale/')
}

# variants info for BIG40
if (!file.exists('data/MR/MR_sumstats/BIG40/variants.txt.gz')&
    !file.exists('data/MR/MR_sumstats/BIG40/variants.txt')){
  system('wget https://open.win.ox.ac.uk/ukbiobank/big40/release2/variants.txt.gz')
  system('mv variants.txt.gz data/MR/MR_sumstats/BIG40/')
}


# Downloaded sumstats -----------------------------------------------------

ls.downloaded = list.files(path = d.output,pattern = '.txt.gz|.tsv.bgz|.tsv.gz',
                                  full.names = T,recursive = T) %>% 
  c(.,gsub('tsv.gz','tsv.bgz',.)) %>% 
  unique

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
  .[!grepl('irnt',.$`Phenotype Code`)]
ls.rm.multiple = ls.single.neale$`Phenotype Code` %>% 
  strsplit(.,split = '_') %>% lapply(FUN=function(x) x[1]) %>% 
  unlist %>% 
  .[duplicated(.)] %>% unique %>% 
  paste0('^',.,'_',collapse = '|')
ls.single.neale = ls.single.neale %>% 
  .[!grepl(ls.rm.multiple,.$`Phenotype Code`),] %>% 
  mutate(field_used = strsplit(.$`Phenotype Code`,split = '_') %>% lapply(FUN=function(x) x[1]) %>% 
           unlist) %>% 
  mutate(wget_command=gsub('-O ',paste0('-O ',d.output,'/'),`wget command`)) %>% 
  mutate(file_loc=wget_command %>% strsplit(.,' -O ') %>% lapply(tail,n=1) %>% unlist) %>% 
  .[!.$file_loc %in% ls.downloaded,]

# BIG40
ls.single.big40 = big40.gwasls %>% 
  .[grepl(ls.single.phewas,.$`UKB.ID`),] %>% 
  mutate(wget_command = paste0('curl -o ',d.output,'/IDP_',`UKB.ID`,'.txt.gz -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/',
                              Pheno,'.txt.gz')) %>% 
  mutate(field_used=UKB.ID) %>% 
  mutate(file_loc=wget_command %>% strsplit(.,' ') %>% 
           lapply(head,n=3) %>% lapply(tail,n=1) %>% unlist) %>% 
  .[!.$file_loc %in% ls.downloaded,]

# Download
ls.single.neale$`wget command` %>% 
  c(.,ls.single.big40$wget_command) %>% 
  as.list %>% 
  pblapply(.,system)


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
  .[!grepl('irnt',.$`Phenotype Code`)]
ls.rm.multiple = ls.multiple.neale$`Phenotype Code` %>% 
  strsplit(.,split = '_') %>% lapply(FUN=function(x) x[1]) %>% 
  unlist %>% 
  .[duplicated(.)] %>% unique %>% 
  paste0('^',.,'_',collapse = '|')
ls.multiple.neale = ls.multiple.neale %>% 
  .[!grepl(ls.rm.multiple,.$`Phenotype Code`),] %>% 
  mutate(field_used = strsplit(.$`Phenotype Code`,split = '_') %>% lapply(head,n=1) %>% 
           unlist) %>% 
  mutate(wget_command=gsub('-O ',paste0('-O ',d.output,'/meta/'),`wget command`)) %>% 
  mutate(file_loc=wget_command %>% strsplit(.,' -O ') %>% lapply(tail,n=1) %>% unlist) %>% 
  .[!.$file_loc %in% ls.downloaded,]

# BIG40
ls.multiple.big40 = big40.gwasls %>% 
  .[grepl(ls.multiple.phewas,.$`UKB.ID`),] %>% 
  mutate(wget_command = paste0('curl -o ',d.output,'/meta/IDP_',`UKB.ID`,'.txt.gz -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/',
                               Pheno,'.txt.gz')) %>% 
  mutate(field_used=UKB.ID) %>% 
  mutate(file_loc=wget_command %>% strsplit(.,' ') %>% 
           lapply(head,n=3) %>% lapply(tail,n=1) %>% unlist) %>% 
  .[!.$file_loc %in% ls.downloaded,]

# Download
ls.multiple.neale$`wget command` %>% 
  c(.,ls.multiple.big40$wget_command) %>% 
  as.list %>% 
  pblapply(.,system)

ls.multiple = ls.multiple.neale %>%
  select(wget_command,field_used,file_loc) %>%
  rbind(ls.multiple.big40[,c('wget_command','field_used','file_loc')])

saveRDS(ls.multiple,file=paste0('data/MR/meta.ls.rds'))
