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
f.web = args[2]
d.meta = args[3]
f.Neale = args[4]
f.BIG40 = args[5]
d.output = args[6]

# Neale list: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=178908679
# BIG40 list: https://open.win.ox.ac.uk/ukbiobank/big40/BIG40-IDPs_v4/IDPs.html
# Convert BIG40 to csv: https://www.convertcsv.com/html-table-to-csv.htm

# f.fields = 'data/MR/MR_pheno_ls_UKB.RData'
# f.web = 'data/MR/meta.ls.rds'
# d.meta = 'data/MR/MR_sumstats/meta'
# f.Neale = 'data/MR/UKBB_GWAS_Imputed_v3_201807.tsv'
# f.BIG40 = 'data/MR/BIG40_GWAS.csv'
# d.output = 'data/MR/MR_sumstats'


# Load data ---------------------------------------------------------------
load(f.fields)

ls.separate.field = readRDS(f.web) %>% 
  mutate(file_loc = gsub('.bgz','.gz',file_loc)) %>% 
  select(-wget_command) 

fields.meta = fields.multiple %>% 
  select(field_tag,field_used_1,field_used_2) %>% 
  left_join(ls.separate.field,by=c('field_used_1'='field_used')) %>% 
  rename(file_1=file_loc) %>% 
  left_join(ls.separate.field,by=c('field_used_2'='field_used')) %>% 
  rename(file_2=file_loc) %>% 
  filter(!is.na(file_1)) %>% 
  mutate(file_1=gsub('.bgz','.gz',file_1),file_2=gsub('.bgz','.gz',file_2))

# Neale gwas
neale.gwasls = fread(f.Neale)
ref.neale = fread(paste0(d.output,'/Neale/variants.tsv.gz')) %>% 
  select(variant,A1=ref,A2=alt,SNP=rsid,info,Freq=AF) %>% 
  filter(info>0.8)

# BIG40 IDP gwas
big40.gwasls = read.csv(f.BIG40,skip = 1,stringsAsFactors = F) %>% 
  select(-X,-X1) %>% 
  mutate(Pheno=sprintf("%04d", Pheno))
ref.big40 = fread(paste0(d.output,'/BIG40/variants.txt.gz')) 


# Load reformatted list ---------------------------------------------------

ls.reformatted = list.files(d.meta,pattern = '_formetal.txt.gz',full.names = T,recursive = 
                              T) %>% 
  c(gsub('_formetal.txt.gz','.txt.gz',.),gsub('_formetal.txt.gz','.tsv.gz',.)) 

# Reformat sumstats for metal ---------------------------------------------

ls.f.IDP = c(fields.meta$file_1,fields.meta$file_2) %>% 
  .[grepl('IDP',.)] %>% unique

ls.f.neale = c(fields.meta$file_1,fields.meta$file_2) %>%
  .[! . %in% ls.f.IDP]

# BIG40 IDP
reformat_idp <- function(tmp.name,tmp.ref,af.info){
  ls.running = read_tsv('data/MR/running_list_reformat',col_names = F) %>% 
    .$X1 %>% as.vector %>% as.character
  if(!tmp.name %in% ls.running){
  write_tsv(data.frame(tmp.name),path='data/MR/running_list_reformat',append = T)
  target.pheno = tmp.name %>% strsplit(.,'/') %>% unlist %>% tail(n=1) %>% 
    gsub('IDP_','',.) %>% 
    gsub('.txt.gz','',.,fixed=T)
  target.fname = tmp.name %>% gsub('.txt.gz','_formetal.txt.gz',.,fixed = T)
  tmp.n = tmp.ref %>% 
    .[.$UKB.ID==target.pheno,'N.all.']
  x = fread(tmp.name)
  x = x %>% 
    mutate(pval=10^(-`pval(-log10)`),n=tmp.n) %>%
    left_join(.,af.info,by=c('chr'='chr','pos'='pos','rsid'='rsid','a1'='a1','a2'='a2')) %>% 
    select(snpid=rsid,chr,bpos=pos,a1=a2,a2=a1,
                   beta,se,pval,freq=af,info,n) %>% 
    filter(info>0.8,chr!='0X') %>% 
    mutate(chr=as.numeric(chr)) 
  write_delim(x,path = target.fname)
  }
  gc()
}

system('touch data/MR/running_list_reformat')
ls.f.IDP %>% .[! . %in% ls.reformatted] %>% as.list %>% 
  pblapply(reformat_idp,tmp.ref=big40.gwasls,af.info=ref.big40) 

# Neale lab gwas
reformat_neale <- function(tmp.name,tmp.ref){
  target.fname = tmp.name %>% gsub('.tsv.gz','_formetal.txt.gz',.,fixed = T)
  x = fread(tmp.name)
  x = x %>% 
    left_join(.,tmp.ref,by='variant') %>% 
    filter(chr!='0X') %>%
    mutate(chr=as.numeric(chr)) %>% 
    select(chr,bpos=pos,snpid=rsid,a1=A1,a2=A2,
           beta,se,z,pval,freq=Freq,n=n_complete_samples) 
      
  write_delim(x,path = target.fname)
}

ls.f.neale %>% .[! . %in% ls.reformatted] %>% as.list %>% 
  pblapply(reformat_neale,tmp.ref=ref.neale) 