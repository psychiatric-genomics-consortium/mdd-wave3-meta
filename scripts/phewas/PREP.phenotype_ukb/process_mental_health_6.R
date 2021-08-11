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
f.mhq = args[3]                     # f.mhq = 'data/2021-04-phenotypes-ukb44797/derived/MHQ.1907.ukb24262.Process_MH_Questionnaire_Output.rds'
f.mdd = args[4]                     # f.mdd = 'data/2021-04-phenotypes-ukb44797/derived/ukb24262-mdd.mdd_phenotypes.rds'
f.touchscreen = args[5]             # f.touchscreen = 'data/2021-04-phenotypes-ukb44797/Touchscreen.rds'
f.output_data = args[6]
f.output_dictionary = args[7]


# Load data ---------------------------------------------------------------

fields.all=fread(f.dic,header=T,stringsAsFactors=F)
ref.coding=read.csv(f.coding,header=T,stringsAsFactors=F) 

touchscreen = readRDS(f.touchscreen)
mdd = readRDS(f.mdd)
mhq = readRDS(f.mhq)

# Touchscreen -------------------------------------------------------------
# Derived data from touchscreen, categories 100060 and 100061 included as data
# Touchscreen data was then merged with online MHQ

ls.fields.available = colnames(touchscreen) %>% .[!. %in% 'f.eid'] %>% 
  strsplit(., '\\.') %>% 
  bind_cols %>% 
  t %>% data.frame(.,stringsAsFactors=F) %>% 
  .$X2 %>% as.numeric

fields.mh.touchscreen = fields.all %>% 
  .[grep('Data',.$ItemType),] %>% 
  .[grep('100061|100060',.$Category),] %>% 
  .[.$FieldID %in% ls.fields.available,] %>% 
  .[!grepl('(pilot)',.$Field),]

# items to remove as covered by total scores 
# (e.g. individual items removed for calculating neuroticism)
ls.fields.rm = seq(1920,2030,10) %>% 
  c(.,seq(2050,2110,10)) %>% 
  c(.,20122:20127)
fields.mh.touchscreen = fields.mh.touchscreen %>% 
  .[!.$FieldID %in% ls.fields.rm,]

fields.mh.touchscreen.block = fields.mh.touchscreen %>% 
  .[grepl('Categorical multiple',.$ValueType),] %>% 
  mutate(category='Mental health',
         field_tag=paste0('f.',FieldID),
         field_used=as.character(FieldID))
fields.mh.touchscreen.nonblock = fields.mh.touchscreen %>% 
  .[!grepl('Categorical multiple',.$ValueType),] %>% 
  mutate(category='Mental health',
         field_tag=paste0('f.',FieldID),
         field_used=as.character(FieldID))

mh.touchscreen = touchscreen %>% 
  select(f.eid, matches(paste0('f.',fields.mh.touchscreen$FieldID,'\\.',collapse = '|'))) 


# Block data: count multiple answers 
# (number of recent negative events, social support and maniac symptoms)

extract_blockpheno <- function(field_dat,dat.dic,instances=0:2){
  FieldID=field_dat$FieldID
  master.dat=get(field_dat$dat)
  
  for (i in instances){
    tmp.block.dat<- master.dat %>%
      select(matches(paste0('f.',FieldID,'\\.'))) %>%
      select(matches(paste0('\\.',instances,'\\.', collapse="|"))) 
    tmp.interpolate.dat=rowSums(tmp.block.dat!='None of the above',na.rm=T)
    tmp.interpolate.dat[rowSums(!is.na(tmp.block.dat))==0]=NA
    if (i==instances[1]){
      tmp.output.dat=data.frame(tmp.interpolate.dat,stringsAsFactors = F)
      colnames(tmp.output.dat)=paste0('f.',FieldID)
    }else{
      loc.to.interpolate=is.na(tmp.output.dat[,names(tmp.output.dat)])
      tmp.output.dat[,names(tmp.output.dat)][loc.to.interpolate]=
        tmp.interpolate.dat[loc.to.interpolate]             
    }
    
    # keep a record of N derived from interpolated data
    if(i==instances[1]){
      count.n.record=sum(!is.na(tmp.interpolate.dat))
    }else{
      count.n.instance=sum(!is.na(tmp.interpolate.dat[loc.to.interpolate]))
      count.n.record=c(count.n.record,count.n.instance)
      }
  } 
  
  tmp.output.dat$f.eid=master.dat$f.eid
  output.pheno=list('phenotype'=tmp.output.dat,'count'=count.n.record)
  return(output.pheno)  
}

field.input=data.frame(FieldID=fields.mh.touchscreen.block$FieldID,
                       dat='mh.touchscreen',stringsAsFactors=F) %>%
  split(., seq(nrow(.)))
mh.touchscreen.blockrmed.dat = pblapply(field.input,FUN=extract_blockpheno,
                                    dat.dic=fields.mh.touchscreen) %>% 
  pblapply(.,function(x) x$phenotype) %>% 
  Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="f.eid"), .)

# Non-block fields: combine multiple instances
extract_pheno <- function(field_dat,dat.dic,instances=0:2){
  FieldID=field_dat$FieldID
  master.dat=get(field_dat$dat)
  tmp.dic=filter(dat.dic,FieldID==FieldID)
  # get data **
  tmp.block.dat <- master.dat %>%
    select(matches(paste0('f.',FieldID,'\\.'))) %>%
    select(matches(paste0('\\.',instances,'\\.', collapse="|")))
  
  # No multiple rounds
  tmp.dat=data.frame(tmp.block.dat[,1])
  colnames(tmp.dat)=paste0('f.',FieldID)
  count.n.record<-data.frame(c(NA,NA,NA))
  rownames(count.n.record)=paste0('instance.',instances)
  colnames(count.n.record)=paste0('f.',FieldID)
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
field.input=data.frame(FieldID=fields.mh.touchscreen.nonblock$FieldID,
                       dat='mh.touchscreen',stringsAsFactors=F) %>%
  split(., seq(nrow(.)))
mh.touchscreen.singleinstance = 
  pblapply(field.input,FUN=extract_pheno,dat.dic=fields.mh.touchscreen)
mh.touchscreen.singleinstance.dat = mh.touchscreen.singleinstance %>% 
  pblapply(.,function(x) x$phenotype) %>% 
  Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="f.eid"), .)

# Numerise ordinal vars
ls.numerise = fields.mh.touchscreen %>% 
  .[.$ValueType %in% c('Categorical single'),] %>% 
  .$FieldID %>% 
  paste0('f.',.)

targetdata = mh.touchscreen.singleinstance.dat
for (i in 2:ncol(targetdata)){
  tmp.dat = targetdata[,i]
  tmp.dat[tmp.dat=='Do no know']=NA
  tmp.dat[tmp.dat=='Prefer not to answer']=NA
  if ((colnames(targetdata)[i] %in% ls.numerise) & is.factor(targetdata[,i])){
    tmp.dat=as.numeric(tmp.dat)
    tmp.dat[tmp.dat<0]=NA
    targetdata[,i]=tmp.dat
  }else if((colnames(targetdata)[i] %in% ls.numerise) & !is.factor(targetdata[,i])){
    tmp.dat[tmp.dat<0]=NA
    targetdata[,i]=tmp.dat
  }else{
    targetdata[,i]=tmp.dat
  }
}

mh.touchscreen.singleinstance.numerised = targetdata

mh.touchscreen = merge(mh.touchscreen.singleinstance.numerised,mh.touchscreen.blockrmed.dat,
                       by='f.eid',all=T)


# MDD definition ----------------------------------------------------------

mdd = mdd %>% 
  rename(MDD_nerves = mdd_nerves, MDD_smith = mdd_smith, MDD_ICD = mdd_icd)

fields.mdd = fields.mh.touchscreen[1:3,] %>% 
  mutate(Field=c('MDD_nerves','MDD_smith','MDD_ICD')) %>%
  mutate(FieldID=Field,Link='https://doi.org/10.1192/bjo.2019.100') %>% 
  mutate_at(vars(-Field,-FieldID,-Link),function(x) x=NA) %>% 
  mutate(Notes='Derived MDD definitions',ItemType='Data') %>% 
  mutate(Participants=colSums(!is.na(mdd[,2:ncol(mdd)]))) %>% 
  mutate(category='Mental health',field_tag=FieldID,field_used=NA)


# MHQ ---------------------------------------------------------------------

mhq = mhq %>% 
  rename(MDD_CIDI=Depressed.Ever)
# remove variables that duplicate with touchscreen items
ls.nonMH = colnames(mhq) %>% 
  grep('^SR',.) %>% 
  min %>% 
  (function(x) x-1) %>% 
  2:. %>% 
  colnames(mhq)[.]

ls.col.keep = colnames(mhq) %>% 
  .[!. %in% ls.nonMH] %>% 
  .[!grepl('No.Info|.Screen|.Response|InterviewDepression',.)] %>% 
  .[!grepl('NoSRConditions|MultipleSRConditions|Comorbidity|MH.Questionnaires',.)] 

mhq = mhq %>% 
  .[,ls.col.keep] %>% 
  data.frame %>% 
  select(f.eid,MDD_CIDI,everything())

fields.mhq = fields.mh.touchscreen[1:(length(ls.col.keep)-1),] %>% 
  mutate(Field=ls.col.keep %>% .[!.%in%'f.eid']) %>%
  mutate(FieldID=Field,Link='https://doi.org/10.1192/bjo.2019.100',ItemType='Data',
      Notes='All phenotypes derived from online mental health questionnaire, using the codes provided by Davis et al. See column Link for the protocol paper') %>% 
  mutate_at(vars(-Field,-FieldID,-Link,-Notes,-ItemType),function(x) x=NA) %>% 
  mutate(Field=gsub("([[:lower:]])([[:upper:]])","\\1 \\2",Field)) %>%
  mutate(Field=gsub('^SR','Self Reported ',Field)) %>% 
  mutate(Field=gsub('\\.',' ',Field)) %>%
  mutate(Participants=colSums(!is.na(mhq[,2:ncol(mhq)]))) %>% 
  mutate(category='Mental health',field_tag=FieldID,field_used=NA)


# Merge data --------------------------------------------------------------

mh.output = mh.touchscreen %>% left_join(.,mdd,by='f.eid') %>% left_join(.,mhq,by='f.eid')

fields.mh.output = rbind(fields.mh.touchscreen.nonblock,fields.mh.touchscreen.block,
                         fields.mdd,fields.mhq) 

update_tag <- function(x,tmp.tag){
  if(sum(grepl('f\\.',tmp.tag))>=1){
    real.tag = colnames(x) %>%
      .[grep(paste0('^',tmp.tag,'\\.|^',tmp.tag,'$'),.)]
  }else{
    real.tag = colnames(x) %>%
      .[grep(paste0('^',tmp.tag,'$'),.)]
  }
  return(real.tag)
}

tag.update = fields.mh.output$field_tag %>%
  as.list %>%
  lapply(.,FUN=update_tag,x=mh.output) %>%
  unlist %>%
  as.character

fields.mh.output$field_tag = tag.update


# Remove columns var=0 ----------------------------------------------------

fields.non0var = mh.output %>% sapply(.,as.numeric) %>%
  apply(.,MARGIN=2,var,na.rm=T) %>% as.vector %>%
  (function(x) x!=0) %>% colnames(mh.output)[.] %>% .[!is.na(.)]
mh.output=mh.output[,fields.non0var]

fields.mh.output = fields.mh.output %>% 
  .[.$field_tag %in% colnames(mh.output),]

# Save data and dictionary ------------------------------------------------

saveRDS(mh.output,file=f.output_data)
write.table(fields.mh.output,
            file=f.output_dictionary,sep='\t',quote=T,row.names=F,col.names=T)

