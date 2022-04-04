# Some gwas managed to finish in one go. Sumstats were stored in gz format data for an entire category.
# These files are too big to read in R. 
# This script is used to separate the sumstats by phenotype using bash.
# Outputs: pheno_<category name>_step2_<pheno name>.regenie
# Bash scripts will be store where the sumstats base


library(dplyr)
library(readr)
library(pbapply)

# Inputs ------------------------------------------------------------------

ls.regenie.gz = list.files('/exports/eddie/scratch/xshen33/phewas_gwas',pattern = '.regenie.gz') %>% 
  .[grep('pheno_',.)]
ls.category = ls.regenie.gz %>% 
  gsub('_step2.regenie.gz','',.) %>% 
  gsub('pheno_','',.)

ls.f.pheno = paste0('/exports/eddie/scratch/xshen33/phewas_gwas/pheno_',ls.category,'.tsv')

dic.pheno_cate = data.frame(cate = ls.category,
                            f.regenie = ls.regenie.gz,
                            f.pheno = ls.f.pheno)

# Function ----------------------------------------------------------------

pheno_col <- function(y.pheno_name,y.pheno_dic,y.data_col){
  y.seq = y.pheno_dic %>% 
    filter(field_id ==y.pheno_name) %>% 
    {.$seq_pheno}
  out.colnames = c('CHROM','GENPOS','ID','ALLELE0','ALLELE1','A1FREQ','INFO','N','TEST',
                   paste0(c('BETA','SE','CHISQ','LOG10P'),'.Y',y.seq),'EXTRA')
  out.colno = which(y.data_col %in% out.colnames)
  return(out.colno)
}

create_bash_cmmd <- function(z.cols){
  z.cmmd = z.cols %>% 
    paste0('$',.,collapse = ', ') %>% 
    paste0(' awk \'{print ',.,'}\'')
  return(z.cmmd)
}

create_inputs <- function(x.category){
  
  # load pheno dat
  x.pheno = dic.pheno_cate %>% filter(cate==x.category) %>% {.$f.pheno} %>% 
    read_tsv(.) %>% 
    select(-FID,-IID)
  
  # create a list of phenos to produce outputs
  x.ls.pheno = colnames(x.pheno)
  x.dic.seq_pheno = data.frame(seq_pheno = 1:length(x.ls.pheno),
                               field_id = x.ls.pheno)
  # target file colunm names
  pheno_s.names = c('BETA','SE','CHISQ','LOG10P') %>% 
    expand.grid(.,x.dic.seq_pheno$seq_pheno) %>% as.data.frame
  pheno_s.names = paste0(pheno_s.names[,1],'.Y',pheno_s.names[,2])
  dic.colnames = c('CHROM','GENPOS','ID','ALLELE0','ALLELE1','A1FREQ','INFO','N','TEST') %>% 
    c(.,pheno_s.names,'EXTRA')
  
  # col no.s for each pheno
  col.no_s = x.dic.seq_pheno$field_id %>% as.list %>% 
    pblapply(pheno_col,y.pheno_dic=x.dic.seq_pheno,y.data_col=dic.colnames) 
  
  # create bash script
  x.regeniefile = dic.pheno_cate %>% filter(cate==x.category) %>% {.$f.regenie} %>% 
    gsub('(','\\(',.,fixed=T) %>% gsub(')','\\)',.,fixed=T)
  x.pheno_outfiles = x.regeniefile %>% gsub('.regenie.gz','_',.) %>% 
    paste0(.,x.dic.seq_pheno$field_id,'.regenie')
  bash_cmmd = col.no_s %>% pblapply(create_bash_cmmd) %>% 
    lapply(FUN=function(m){
      paste0('zcat ',x.regeniefile,' |',m,' > ')
    }) %>% unlist %>% 
    paste0(.,x.pheno_outfiles)
  return(bash_cmmd)
  # write script
  
}


# Generate bash scripts ---------------------------------------------------

bash_script = ls.category %>% as.list %>% pblapply(create_inputs) %>% unlist
bash_script = data.frame(cmmd = bash_script)

write_tsv(bash_script,file='/exports/eddie/scratch/xshen33/phewas_gwas/reformat_singlerun.sh',col_names = F)