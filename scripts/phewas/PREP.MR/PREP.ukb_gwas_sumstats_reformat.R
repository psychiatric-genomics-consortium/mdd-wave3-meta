library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

d.dat = args[1]
d.output = args[2]

# d.dat = '/exports/eddie/scratch/xshen33/phewas_gwas'
# d.output = 'data/MR/MR_sumstats'


# Load phenos -------------------------------------------------------------

ls.category = list.files(path=d.dat,pattern='.regenie.gz') %>% 
  .[grepl('step2',.)] %>% 
  gsub('.regenie.gz','',.,fixed = T) %>% 
  gsub('_([0-9])','',.) %>% 
  unique

system('touch data/MR/running_regenie_reformat')

get_pheno <- function(pheno.name,tmp.dat,ref.pheno){
  real_field = ref.pheno[pheno.name,'field_name']
  new.fname = paste0(real_field,'.regenie.gz')
  
  tmp.new.dat = tmp.dat %>% 
    select(-matches('.Y'),matches(pheno.name)) %>% 
    rename_with(~gsub(paste0('.',pheno.name), "", .x, fixed = TRUE)) %>% 
    mutate(p=10^(-LOG10P),Freq1=1-A1FREQ) %>% 
    select(SNP=ID,Allele1=ALLELE1,Allele2=ALLELE0,Effect=BETA,se=SE,
           p,Freq1,CHR=CHROM,BP=GENPOS,N) %>% 
    as.data.frame
  
  system(paste0('echo ',new.fname,'>> data/MR/running_regenie_reformat'))
  write_tsv(tmp.new.dat,paste0(d.output,'/',tmp.field,'.regenie.gz'))
  
}

load_sumstats <- function(tmp.field,tmp.path){
  ls.processed = list.files(d.output,pattern='regenie.gz') %>% 
    gsub('.regenie.gz','',.,fixed = T)
  ls.running = read_tsv('data/MR/running_regenie_reformat',col_names=F) %>% 
    .$X1 %>% as.character %>% as.vector
  ls.avoid = c(ls.processed,ls.running)
  
  if(sum(tmp.field %in% ls.avoid)==0){
    
    ls.pheno = list.files(path=tmp.path,pattern=tmp.field,full.names = T) %>% 
      .[grepl('step2',.)] %>% 
      .[grepl('regenie.Ydict',.)] %>% 
      read_delim(.,col_names=F) %>% as.data.frame %>% 
      rename(pheony_field=X1,field_name=X2)
    rownames(ls.pheno)=ls.pheno$pheony_field
    
    tmp.fs = list.files(path=tmp.path,pattern=tmp.field,full.names = T) %>% 
      .[grepl('step2',.)] %>% 
      .[grepl('regenie.gz',.)]
    
    dat.gwas = tmp.fs %>% as.list %>% 
      lapply(read_delim,delim=" ") %>% 
      bind_rows %>% 
      as.data.frame
    
    ls.pheno$pheony_field %>% 
      as.list %>% 
      lapply(get_pheno,tmp.dat=dat.gwas,ref.pheno)

  }
  gc()
}


ls.category %>% as.list %>% 
  pblapply(load_sumstats,tmp.path=d.dat)
