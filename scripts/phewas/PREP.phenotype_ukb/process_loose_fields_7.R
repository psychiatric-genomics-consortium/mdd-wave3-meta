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
f.pheno_dir = args[3]               # f.pheno_dir = 'data/2021-04-phenotypes-ukb44797/Imaging.rds'
f.categories = args[4]              # f.categories = 'results/phewas/data_dictionary/category_loose_fields'
f.output_data = args[5]
f.output_dictionary = args[6]



# Load inputs -------------------------------------------------------------

fields.loose=fread(f.loose_field,header=T,stringsAsFactors=F) 
# Correct a typo in the data dictionary
fields.loose[fields.loose$FieldID==23342,'Field']='Femur wards bone area (right)' 
ref.coding=read.csv(f.coding,header=T,stringsAsFactors=F) 
ref.category=fread(f.categories,header=F,stringsAsFactors=F) %>% 
  rename(Category=V1,category=V2) %>% data.frame

g.path = f.pheno_dir %>% strsplit(.,'/') %>% unlist %>% 
  .[1:(length(.)-1)] %>% 
  paste0(.,collapse = '/')


# Load data ---------------------------------------------------------------

# Load available data fields from files ***
fs = list.files(path=g.path,full.names = T)%>% 
  (function(x) x[grep('\\.rds',x)])

fs.list=fs %>% as.list
# extracting available fields in the files
extract_available_fields <- function(fname,dat.dic,path=g.path){
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

# fields.loose %>% 
#   .[as.numeric(.$FieldID) %in% as.numeric(fields.available$field),] %>% 
#   write.csv(.,file='data/tmp_available_loose_field',quote=T,row.names = F,col.names = T)

for(x in files.to.obj){
  tmp.dat=readRDS(x$file)
  tmp.fields=c('f.eid',paste0('f.',fields.available$field,'\\.'))
  tmp.fields=paste0(tmp.fields,collapse='|')
  tmp.dat<-tmp.dat %>%
    select(matches(tmp.fields))
  
  eval(parse(text=paste0(x$obj.name,'=tmp.dat')))
  cat(paste0(x$obj.name,'\n'))
}

# Merge data
loose_field.dat = ls(pattern = '.rds') %>% as.list %>% lapply(get) %>% 
  Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="f.eid"), .)

# Block data --------------------------------------------------------------

fields.loose.block = fields.loose %>% 
  .[.$FieldID %in% as.numeric(fields.available$field),] %>% 
  .[grepl('Categorical multiple',.$ValueType),] %>% 
  .[grepl('problems|bone site(s)|diagnosed|Medication|medication|supplements|Types|Location of pain',
          .$Field),] %>% 
  mutate(field_tag=paste0('f.',FieldID),field_used=as.character(FieldID))

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

field.input=data.frame(FieldID=fields.loose.block$FieldID,
                       dat='loose_field.dat',stringsAsFactors=F) %>%
  split(., seq(nrow(.)))

loose_field.blockrmed.dat = pblapply(field.input,FUN=extract_blockpheno,
                                        dat.dic=fields.loose.block) %>% 
  pblapply(.,function(x) x$phenotype) %>% 
  Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="f.eid"), .)


# Non-block data ----------------------------------------------------------

fields.loose.nonblock = fields.loose %>% 
  .[.$FieldID %in% as.numeric(fields.available$field),] %>%
  .[!.$FieldID %in% fields.loose.block$FieldID,] 

extract_pheno <- function(field_dat,dat.dic,instances=0:2){
  FieldID=field_dat$FieldID
  master.dat=get(field_dat$dat)
  tmp.dic=filter(dat.dic,FieldID==FieldID)
  # get data **
  tmp.block.dat <- master.dat %>%
    select(matches(paste0('f.',FieldID,'\\.'))) %>%
    select(matches(paste0('\\.',instances,'\\.', collapse="|"))) 
  
  # No multiple rounds
  tmp.dat=data.frame(tmp.block.dat[,1],stringsAsFactors = F)
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
field.input=data.frame(FieldID=fields.loose.nonblock$FieldID,
                       dat='loose_field.dat',stringsAsFactors=F) %>%
  split(., seq(nrow(.)))
loose_field.singleinstance = 
  pblapply(field.input,FUN=extract_pheno,dat.dic=fields.loose.nonblock)
loose_field.singleinstance.dat = loose_field.singleinstance %>% 
  pblapply(.,function(x) x$phenotype) %>% 
  Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="f.eid"), .)

# Numerise ordinal vars
ls.numerise = fields.loose.nonblock %>% 
  .[.$ValueType %in% c('Categorical single'),] %>% 
  .$FieldID %>% 
  paste0('f.',.)

targetdata = loose_field.singleinstance.dat
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

loose_field.singleinstance.numerised = targetdata 


# Combine left n right data -----------------------------------------------

find_right<-function(x,tmp.ref){
  left.Field=x$Field
  right.Field=gsub('left','right',left.Field)
  left.Category=x$Category
  right.ref=tmp.ref %>% 
    filter(Category==left.Category) %>% 
    filter(Field==right.Field)
  if(nrow(right.ref)>1){cat('Too many matching fields')}else if(nrow(right.ref)==1){
    return(right.ref)
  }
}

left.ref=fields.loose.nonblock[grep('left',fields.loose.nonblock$Field),]

right.ref = left.ref %>% 
  split(.,seq(nrow(.))) %>% 
  lapply(.,find_right,tmp.ref=fields.loose.nonblock) %>% 
  bind_rows

left.fields=left.ref$FieldID
right.fields=right.ref$FieldID


targetdata=loose_field.singleinstance.numerised
left.dat=targetdata[,paste0('f.',left.fields)]
right.dat=targetdata[,paste0('f.',right.fields)]

for (i in 1:ncol(left.dat)){
  aveg.source=cbind(left.dat[,i],right.dat[,i])
  tmp.aveg.dat=data.frame(f.eid=targetdata$f.eid,rowMeans(aveg.source,na.rm=T))
  tmp.aveg.dat[rowSums(!is.na(aveg.source))==0,2]=NA
  colnames(tmp.aveg.dat)[2]=paste0(colnames(left.dat)[i],'.lnr')
  if (i==1){
    aveg.dat=tmp.aveg.dat
  }else{
    aveg.dat=merge(aveg.dat,tmp.aveg.dat,by='f.eid',all.x=T)
  }
}

rm.fields=c(right.fields,left.fields)
rm.fields=paste0('f.',rm.fields)

loose_field.singleinstance.numerised.lnr=loose_field.singleinstance.numerised %>% 
  .[,!colnames(.) %in% rm.fields] %>% 
  left_join(.,aveg.dat,by='f.eid')

fields.loose.nonblock = fields.loose.nonblock %>% 
  .[!.$FieldID %in% right.fields,] %>% 
  mutate(field_tag=paste0('f.',FieldID),field_used=as.character(FieldID))

lnr.replacement.input = data.frame(FieldID=left.fields,
                             used=paste0(left.fields,',',right.fields),
                             stringsAsFactors=F) 

fields.loose.nonblock.lnr = fields.loose.nonblock 
fields.loose.nonblock.lnr[match(left.fields,fields.loose.nonblock.lnr$FieldID),'field_used']=
  paste0(left.fields,',',right.fields)



# Merge data --------------------------------------------------------------

loose_field.output = loose_field.blockrmed.dat %>% 
  left_join(.,loose_field.singleinstance.numerised.lnr,by='f.eid') 

# New data dictionary
fields.loose.output = rbind(fields.loose.block,fields.loose.nonblock.lnr) %>%
  left_join(.,ref.category,by='Category') %>% 
  select(Path,Category,everything(),category,field_tag,field_used)

update_tag <- function(x,tmp.tag){
  real.tag = colnames(x) %>%
    .[grep(paste0('^',tmp.tag,'\\.|^',tmp.tag,'$'),.)]
  return(real.tag)
}

tag.update = fields.loose.output$field_tag %>%
  as.list %>%
  lapply(.,FUN=update_tag,x=loose_field.output) %>%
  unlist %>%
  as.character

fields.loose.output$field_tag = tag.update


# Remove columns var=0 ----------------------------------------------------


fields.non0var = loose_field.output %>% sapply(.,as.numeric) %>%
  apply(.,MARGIN=2,var,na.rm=T) %>% as.vector %>%
  (function(x) x!=0) %>% colnames(loose_field.output)[.] %>% .[!is.na(.)]
loose_field.output=loose_field.output[,fields.non0var]

fields.loose.output = fields.loose.output %>% 
  .[.$field_tag %in% colnames(loose_field.output),]


# Save data and dictionary ------------------------------------------------

saveRDS(loose_field.output,file=f.output_data)
write.table(fields.loose.output,
            file=f.output_dictionary,sep='\t',quote=F,row.names=F,col.names=T)

