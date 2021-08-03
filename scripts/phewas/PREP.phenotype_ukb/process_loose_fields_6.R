library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.loose_field = args[1]             # f.loose_field = 'data/loose_fields'
f.coding = args[2]                  # f.coding = 'data/Codings.csv'
f.pheno_dir = args[3]               # f.pheno_dir = 'data/2021-04-phenotypes-ukb44797'
f.output_data = args[4]
f.output_dictionary = args[5]



# Load inputs -------------------------------------------------------------

fields.loose=fread(f.loose_field,header=T,stringsAsFactors=F)
ref.coding=read.csv(f.coding,header=T,stringsAsFactors=F) 

g.path = f.pheno_dir


# Load data ---------------------------------------------------------------

# Load available data fields from files ***
fs = list.files(path=g.path,recursive=T,full.names = T)%>% 
  (function(x) x[grep('\\.rds',x)])

fs.list=fs %>% as.list
# extracting available fields in the files
extract_available_fields <- function(fname,dat.dic,path='/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/'){
  tmp.dat=readRDS(fname)
  ls.field=colnames(tmp.dat)[2:ncol(tmp.dat)] %>%
    strsplit(.,'\\.') %>%
    lapply(function(x) x[2]) %>%
    unlist %>%
    unique
  keep.field=ls.field[ls.field %in% dat.dic$FieldID]
  if (length(keep.field>0)){
    ls.field.output=data.frame(field=keep.field,file=fname,stringsAsFactors=F)    
  }else{
    ls.field.output=data.frame(field=NA,file=fname,stringsAsFactors=F)
  }
  
  return(ls.field.output)
}      

fields.available=pblapply(fs.list,extract_available_fields,dat.dic=fields.loose) %>% 
  bind_rows %>% filter(.,!is.na(field)) 
fields.available$obj.name<- fields.available$file %>%
  strsplit(.,'/') %>%
  lapply(function(x) x[length(x)]) %>%
  unlist

# Load data needed ***
files.to.obj=filter(fields.available,!duplicated(file))%>%
  split(., seq(nrow(.)))

for(x in files.to.obj){
  tmp.dat=readRDS(x$file)
  tmp.fields=c('f.eid',paste0('f.',fields.available$field,'\\.'))
  tmp.fields=paste0(tmp.fields,collapse='|')
  tmp.dat<-tmp.dat %>%
    select(matches(tmp.fields))
  
  eval(parse(text=paste0(x$obj.name,'=tmp.dat')))
  cat(paste0(x$obj.name,'\n'))
}


