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

ls.pheno = list.files(path=d.dat,pattern='.regenie$') %>% 
  .[grepl('step2',.)] %>% 
  gsub('.regenie','',.,fixed = T) %>% 
  strsplit(split = '_') %>% 
  lapply(tail,n=1) %>% unlist %>% unique


system('touch data/MR/running_regenie_reformat')

get_pheno <- function(tmp.dat){
  tmp.new.dat = tmp.dat %>% 
    rename_with(~gsub('.Y[0-9]+', '', .x)) %>% 
    mutate(p=10^(-LOG10P),Freq1=1-A1FREQ) %>% 
    select(SNP=ID,Allele1=ALLELE1,Allele2=ALLELE0,Effect=BETA,se=SE,
           p,Freq1,CHR=CHROM,BP=GENPOS,N) %>% 
    as.data.frame
  return(tmp.new.dat)
}

load_sumstats <- function(tmp.field,tmp.path){
  ls.processed = list.files(d.output,pattern='regenie.gz') %>% 
    gsub('.regenie.gz','',.,fixed = T)
  ls.running = read_tsv('data/MR/running_regenie_reformat',col_names=F) %>% 
    .$X1 %>% as.character %>% as.vector
  ls.avoid = c(ls.processed,ls.running)
  
  system(paste0('echo ',tmp.field,'>> data/MR/running_regenie_reformat'))
  
  if(sum(tmp.field %in% ls.avoid)==0){
    dat.pheno = list.files(path=tmp.path,pattern=paste0(tmp.field,'.regenie$'),full.names = T) %>% as.list %>% 
      pblapply(read_delim,col_names=T) %>% lapply(get_pheno) %>% 
      bind_rows
    new.fname = paste0(d.output,'/',tmp.field,'.regenie.gz')
    write_tsv(dat.pheno,file=new.fname)
  }
  gc()
}


ls.pheno %>% as.list %>% 
  pblapply(load_sumstats,tmp.path=d.dat)
