library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dat.brain = args[1]
f.dat.cog = args[2]
f.dat.diet = args[3]
f.dat.activity = args[4]
f.dat.mental = args[5]
f.dat.loose = args[6]
f.dat.cov = args[7]

f.dic_dir = args[8]             
f.PRS = args[9]
f.localgwas = args[10]
d.output = args[10]

# 
# f.dat.brain = 'data/dat.imaging_chunk.rds'
# f.dat.cog = 'data/dat.cognition_chunk.rds'
# f.dat.diet = 'data/dat.diet_chunk.rds'
# f.dat.activity = 'data/dat.activity_chunk.rds'
# f.dat.mental = 'data/dat.mental_health_chunk.rds'
# f.dat.loose = 'data/dat.loose_field_chunk.rds'
# f.dat.cov = 'data/dat.addional_covariates_chunk.rds'
# f.dic_dir = 'results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt'
# f.PRS = 'data/PRS_all.rds'
# f.localgwas = 'data/MR/local.gwas.rds'
# d.output = '/exports/eddie/scratch/xshen33/phewas_gwas'


# Load data ---------------------------------------------------------------

fs = ls(pattern = 'f.dat') %>% as.list %>% lapply(get)

phewas.dat = fs %>% lapply(readRDS) %>% 
  Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="f.eid",all=T), .) %>% 
  right_join(.,readRDS(f.PRS),by='f.eid')

age_sex = phewas.dat %>% 
  select(f.eid, cov.sex,cov.age_instance0)

pa_covs = phewas.dat %>% 
  select(f.eid,f.90051.0.0,f.90161.0.0,f.90162.0.0,
         f.90163.0.0,f.90164.0.0,f.90165.0.0,f.90166.0.0,
         f.90167.0.0,f.90168.0.0,f.90169.0.0)

bolt_covars = read_tsv('/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/gwas/BOLT/whitebritish_centre_array_flashpcs_457k.tsv.gz') %>% 
  select(FID:PC8) %>% 
  as.data.frame %>% 
  left_join(.,age_sex,by=c('FID'='f.eid')) %>% 
  sapply(.,as.numeric) %>% 
  as.data.frame

bolt_covars_pa = bolt_covars %>% 
  left_join(.,pa_covs,by=c('FID'='f.eid')) %>% 
  sapply(.,as.numeric) %>% 
  as.data.frame

ls.gwas = readRDS(f.localgwas)



# Extract pheno data ------------------------------------------------------

prep_pheno <- function(x,tmp.category,tmp.ref){
  block.ref = tmp.ref %>% 
    filter(category==tmp.category)
  ls.pheno = c('f.eid',block.ref$field_tag)
  x.save = x %>% 
    .[,ls.pheno]
  
  ls.sub = rowSums(!is.na(x.save))>1
  x.save = x.save %>% .[ls.sub,] %>% 
    mutate(IID=f.eid) %>% 
    rename(FID=f.eid) %>% 
    select(FID,IID,everything()) %>% 
    sapply(as.numeric) %>% 
    as.data.frame
  return(x.save)
}

pheno.input = unique(ls.gwas$category) %>% as.list %>% 
  pblapply(.,prep_pheno,x=phewas.dat,tmp.ref=ls.gwas)

names(pheno.input) = unique(ls.gwas$category) %>% gsub(' ','_',.)

# regenie pipeline ---------------------------------------------------------

write_tsv(bolt_covars, paste0(d.output,'/phewas_gwas_covariates.tsv'))
write_tsv(bolt_covars_pa, paste0(d.output,'/pa_gwas_covariates.tsv'))

names(pheno.input) %>% as.list %>% 
  pblapply(.,FUN=function(x) write_tsv(pheno.input[[x]], paste0(d.output,'/pheno_',x,'.tsv')))

# run gwas: https://rgcgithub.github.io/regenie/recommendations/
# Check regenie scripts

# Create regenie input list -----------------------------------------------

# regenie step 1 (use genotyped data)
inputs.step1 = data.frame(f.pheno = paste0(d.output,'/pheno_',names(pheno.input),'.tsv')) %>% 
  mutate(f.cov = ifelse(grepl('pain',f.pheno),paste0(d.output,'/pa_gwas_covariates.tsv'),
                        paste0(d.output,'/phewas_gwas_covariates.tsv'))) %>% 
  mutate(f.prefix = gsub('.tsv','_tmp',f.pheno)) %>% 
  mutate(f.out = gsub('.tsv','_step1',f.pheno))

# regenie step 2 (use imputed data separated by chromosome)
inputs.step2 = inputs.step1 %>% 
  select(f.pheno,f.cov,f.out) %>% 
  mutate(f.pred = paste0(f.out,'_pred.list')) %>% 
  mutate(f.out = gsub('_step1','_step2',f.out)) %>% 
  select(f.pheno,f.cov,f.pred,f.out)

write_tsv(inputs.step1, 'data/MR/regenie_step1_ls.tsv',col_names = F)
write_tsv(inputs.step2, 'data/MR/regenie_step2_ls.tsv',col_names = F)



# Create regenie input list for interrupted GWAS (step 2) -----------------

# update output name
update_OutputName <- function(output.pattern){
  tmp.pattern = output.pattern %>% basename
  tmp.dir = output.pattern %>% gsub(tmp.pattern,'',.)
  
  ls.log = list.files(path=tmp.dir,pattern=tmp.pattern,full.names = F) %>% .[grep(pattern='.log',.)]
  
  if(length(ls.log)==0){
    final.output = paste0(output.pattern,'_1')
  }else{
    ls.block_log = ls.log %>% gsub(tmp.pattern,'',.) %>% gsub('_','',.) %>% gsub('.log','',.)
    if(length(ls.block_log)==1&nchar(ls.block_log)==0){
      final.n_block=1+1
    }else{
      final.n_block=ls.block_log %>% as.numeric %>% {max(.,na.rm=T)} %>% {.+1}
    }
    final.output = paste0(output.pattern,'_',final.n_block)
  }
  return(final.output)
}

inputs.step2 = inputs.step2 %>% 
  mutate(f.out_old=f.out) %>% 
  mutate(f.out=f.out %>% gsub('_[0-9]','',.) %>% 
           as.list %>% 
           lapply(update_OutputName) %>% unlist)

# remove those GWAS that finished running, identify starting block
find_finishedGWAS <- function(output.pattern){
  tmp.pattern = output.pattern %>% basename
  tmp.dir = output.pattern %>% gsub(tmp.pattern,'',.)
  
  ls.log = list.files(path=tmp.dir,pattern=tmp.pattern,full.names = F) %>% .[grep(pattern='.log',.)]
  if(length(ls.log)==0){
    starting.block = 1
  }else{
    ls.block_log = ls.log %>% gsub(tmp.pattern,'',.) %>% gsub('_','',.) %>% gsub('.log','',.)
    if(length(ls.block_log)==1&nchar(ls.block_log)==0){
      tmp.log = read_tsv(paste0(tmp.dir,'/',ls.log),col_names=F) 
    }else{
      final.n_block=ls.block_log %>% as.numeric %>% {max(.,na.rm = T)}
      tmp.log = read_tsv(paste0(tmp.dir,'/',ls.log[grep(paste0('_',final.n_block),ls.log)]),col_names=F)
    }
    tmp.log = tmp.log$X1 %>% 
    .[grepl('block \\[',.)] %>% 
    .[grepl('done',.)] %>% 
    strsplit(.,split = '\\]') %>% lapply(.,FUN=function(x) head(x,n=1)) %>% unlist %>% 
    strsplit(.,split = '\\[') %>% lapply(.,FUN=function(x) x[2]) %>% unlist %>% tail(n=1) %>% 
    strsplit(.,split = '/') %>% unlist %>% as.numeric
    
    if(tmp.log[1]==tmp.log[2]){starting.block=NA}else{
      starting.block=tmp.log[1]+1
    }
  }
  return(starting.block)
}

inputs.step2 = inputs.step2 %>% 
  mutate(starting.block=inputs.step2$f.out_old %>% gsub('_[0-9]','',.) %>% 
           as.list %>% 
           lapply(find_finishedGWAS) %>% unlist)

inputs.step2.remaining = inputs.step2 %>% 
  filter(!is.na(starting.block)) %>% 
  select(-f.out_old)

write_tsv(inputs.step2.remaining, 'data/MR/regenie_step2_ls_remaining.tsv',col_names = F)

