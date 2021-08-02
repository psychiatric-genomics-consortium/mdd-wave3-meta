library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dic = args[1]                     # f.dic = 'data/Data_Dictionary_Showcase.csv'
f.coding = args[2]                  # f.coding = 'data/Codings.csv'
f.diet = args[3]                    # f.diet = 'data/2021-04-phenotypes-ukb44797/DietByHourRecall.rds'
f.touchscreen = args[4]             # f.touchscreen ='data/2021-04-phenotypes-ukb44797/Touchscreen.rds'
f.output_data = args[5]
f.output_dictionary = args[6]


# Load data ---------------------------------------------------------------

fields.all=read.csv(f.dic,header=T,stringsAsFactors=F)
ref.coding=read.csv(f.coding,header=T,stringsAsFactors=F) 
fields.diet = fields.all %>% 
  .[grep('Online follow-up > Diet|Touchscreen > Lifestyle and environment > Diet|Touchscreen > Lifestyle and environment > Alcohol',.$Path),]

touchscreen=readRDS(f.touchscreen)
diet=readRDS(f.diet) %>% left_join(.,touchscreen,by='f.eid')

ls.fields.available = colnames(diet) %>% .[!. %in% 'f.eid'] %>% 
  strsplit(., '\\.') %>% 
  bind_cols %>% 
  t %>% data.frame(.,stringsAsFactors=F) %>% 
  .$X2 %>% as.numeric

fields.diet = fields.diet %>% 
  .[.$FieldID %in% ls.fields.available,]

# Process -----------------------------------------------------------------

typical_diet_online_label=100020

# Remove fields for pilot data and those that contain texts and types
remove_category <- function(x,kw,col.n,exact=F){
  if (exact==T){
    new.x=x[!grepl(paste0('^',kw,'$'),x[,col.n]),]
  }else{
    new.x=x[!grepl(kw,x[,col.n]),]
  }
  return(new.x)
}

ls.fields.available = colnames(diet) %>% .[!. %in% 'f.eid'] %>% 
  strsplit(., '\\.') %>% 
  bind_cols %>% 
  t %>% data.frame(.,stringsAsFactors=F) %>% 
  .$X2 %>% as.numeric %>% unique
  
fields.diet = fields.diet %>%
  remove_category(.,kw='(pilot)',col.n='Field') %>%
  remove_category(.,kw='Reason',col.n='Field') %>%
  remove_category(.,kw='completed|started|Duration of questionnaire',col.n='Field') %>%
  (function(x) x[,!colnames(x) %in% c(100020,20085,20086)]) %>%
  remove_category(.,kw='type',col.n='Field') %>%
  remove_category(.,kw='Type',col.n='Field') %>% 
  .[.$FieldID %in% ls.fields.available,]

ls.fields.diet = paste0('f.',fields.diet$FieldID,'.',collapse = '|') %>% 
  paste0('f.eid|',.)
diet.output = diet %>% .[,grepl(ls.fields.diet,colnames(.))] %>% 
  .[,grepl('.0$|f.eid',colnames(.))]


# Combine multiple instances
extract_pheno <- function(field_dat,dat.dic,instances=0:2){
  field.id=field_dat$field.id
  master.dat=get(field_dat$dat)
  tmp.dic=filter(dat.dic,FieldID==field.id)
  # get data **
  tmp.block.dat <- master.dat %>%
    select(matches(paste0('f.',field.id,'\\.'))) %>%
    select(matches(paste0('\\.',instances,'\\.', collapse="|")))
  
  # No multiple rounds
  tmp.dat=data.frame(tmp.block.dat[,1])
  colnames(tmp.dat)=paste0('f.',field.id)
  count.n.record<-data.frame(c(NA,NA,NA))
  rownames(count.n.record)=paste0('instance.',instances)
  colnames(count.n.record)=paste0('f.',field.id)
  for (i in 1:ncol(tmp.block.dat)){
    tmp.interpolate.dat=tmp.block.dat[,i]
    loc.to.interpolate=is.na(tmp.dat[,names(tmp.dat)])
    tmp.dat[,names(tmp.dat)][loc.to.interpolate]=tmp.interpolate.dat[loc.to.interpolate]
    # keep a record of N derived from interpolated data
    if (i==1){
      count.n.instance=sum(!is.na(tmp.interpolate.dat))
    }else{count.n.instance=sum(!is.na(tmp.interpolate.dat[loc.to.interpolate]))}
    count.n.record[i,1]=count.n.instance
  } 
  
  tmp.dat$f.eid=master.dat$f.eid
  output.pheno=list('phenotype'=tmp.dat,'count'=count.n.record)
  return(output.pheno)  
}
field.input=data.frame(field.id=fields.diet$FieldID,dat='diet.output',stringsAsFactors=F) %>%
  split(., seq(nrow(.)))
diet.singleinstance = pblapply(field.input,FUN=extract_pheno,dat.dic=fields.diet)

diet.singleinstance.dat = diet.singleinstance %>% 
  pblapply(.,function(x) x$phenotype) %>% 
  Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="f.eid"), .)

# Convert ordinal categorical variables to numeric, rm dubious values for integers

ls.remain.cate = fields.diet %>% 
  .[.$Coding %in% c(100010,100013,100400),] %>% 
  .$FieldID %>% 
  paste0('f.',.)

targetdata = diet.singleinstance.dat
for (i in 2:ncol(targetdata)){
  tmp.dat = targetdata[,i]
  tmp.dat[tmp.dat=='Do no know']=NA
  tmp.dat[tmp.dat=='Prefer not to answer']=NA
  if ((!colnames(targetdata)[i] %in% ls.remain.cate) & is.factor(targetdata[,i])){
    tmp.dat=as.numeric(tmp.dat)
    tmp.dat[tmp.dat<0]=NA
    targetdata[,i]=tmp.dat
  }else if((!colnames(targetdata)[i] %in% ls.remain.cate) & !is.factor(targetdata[,i])){
    tmp.dat[tmp.dat<0]=NA
    targetdata[,i]=tmp.dat
  }else{
    targetdata[,i]=tmp.dat
  }
}

ls.100400 = fields.diet %>% 
  .[.$Coding %in% 100400,] %>% 
  .$FieldID

for (i in ls.100400){
  i = paste0('f.',i)
  tmp.dat=as.character(targetdata[,i])
  tmp.dat[tmp.dat=='Prefer not to answer']=NA
  tmp.dat[tmp.dat=='Yes, because of illness']='Yes'
  tmp.dat[tmp.dat=='Yes, because of other reasons']='Yes'
  tmp.dat=factor(tmp.dat,ordered = T,levels=c('No','Yes'))
  targetdata[,i]=tmp.dat
}

diet.singleinstance.numerised = targetdata


# Output ------------------------------------------------------------------

field.output = fields.diet %>% 
  mutate(field_tag=FieldID,field_used=paste0(FieldID,'.0-2.0'))

saveRDS(diet.singleinstance.numerised,file=f.output_data)
write.table(field.output,file=f.output_dictionary,sep='\t',quote=F,row.names=F,col.names=T)
