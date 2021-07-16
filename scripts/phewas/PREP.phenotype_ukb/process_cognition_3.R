library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.cog_ac = args[1]                 # file location for cognition results (assessment centre)
f.cog_online = args[2]             # file location for cognition results (online)
f.cogvar_QCed = args[3]            # file location for QCed variables (non-outcome vars removed, check Extract field list section)
f.output_data = args[4]
f.output_dictionary = args[5]


# Load data ---------------------------------------------------------------

cognition=readRDS(f.cog_ac)
cognition.2=readRDS(f.cog_online) 
overlap.cols = intersect(colnames(cognition),colnames(cognition.2)) %>% .[!.%in%'f.eid']
cognition.2=cognition.2 %>% 
  .[,! colnames(.) %in% overlap.cols]

cognition=merge(cognition,cognition.2,by='f.eid',all.x=T)

ls.fields.available = colnames(cognition) %>% .[!. %in% 'f.eid'] %>% 
  strsplit(., '\\.') %>% 
  bind_cols %>% 
  t %>% data.frame(.,stringsAsFactors=F) %>% 
  .$X2 %>% as.numeric

fields.cog.outcome=read.delim(f.cogvar_QCed,header=T,stringsAsFactors=T) %>% 
  .[.$FieldID %in% ls.fields.available,]



  

# Extract task list and fields to include  --------------------------------------------------------------------
# remove_category <- function(x,kw,col.n,exact=F){
#   if (exact==T){
#     new.x=x[!grepl(paste0('^',kw,'$'),x[,col.n]),]
#   }else{
#     new.x=x[!grepl(kw,x[,col.n]),]
#   }
#   return(new.x)
# }
# 
# fields.cog = fields.cog %>%
#   remove_category(.,kw='(pilot)',col.n='Field')
# 
# ls.cog.task <- unique(fields.cog$Path) %>% 
#   gsub('UK Biobank Assessment Centre > ','',.) %>%
#   gsub('Online follow-up > ','',.)

# write.table(fields.cog,file='data/data_dictionary/fields.cog_for_process.txt',row.names=F,col.names=T,sep='\t') # for manual check
# remove unnecessary columns (non-outcome variables) manually

# Process block data (multiple rounds) ----------------------------------------------------------------------------

cog_functions_raw = 
  cognition %>%
  select(matches(paste(paste0('f.',fields.cog.outcome$FieldID), collapse="|")))

# find unique list of task * instances
ls.task_instance = colnames(cog_functions_raw) %>% 
  .[ !. %in% 'f.eid'] %>% 
  strsplit(., '\\.') %>% 
  bind_cols %>% 
  t
rownames(ls.task_instance)=NULL
ls.fields_findblock <- apply(ls.task_instance[,1:3],1,paste,collapse=".") %>% unique %>% paste0(.,'.')

# collapse block data
# find fields that have blocks 
# n.fields_block = apply(data.frame(f=ls.fields_findblock),1, function(f) length(grep(f,colnames(cog_functions_raw))))

collapse_block <- function(x,f){
  tmp.block.dat = x[,grep(f,colnames(x))]
  coln = f
  if(is.vector(tmp.block.dat)){
    tmp.outcome = data.frame(tmp.block.dat)
  }else{
    ls.na=rowSums(!is.na(tmp.block.dat))
    tmp.outcome=rowMeans(tmp.block.dat,na.rm=T)
    tmp.outcome[ls.na==0]=NA
    tmp.outcome = data.frame(tmp.outcome)
  }
  colnames(tmp.outcome) = coln
  return(tmp.outcome)
}


cog_functions_noblock <- cog_functions_raw %>%
  sapply(.,as.numeric)  # Don't pipe up here: not enough memory
cog_functions_noblock <- cog_functions_noblock %>% 
  pblapply(as.list(ls.fields_findblock),collapse_block,x=.) %>%
  bind_cols
cog_functions_noblock=data.frame(f.eid=cognition$f.eid,cog_functions_noblock)



# Process data of different instances ----------------------------------------------------------------------------

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
field.input=data.frame(field.id=fields.cog.outcome$FieldID,dat='cog_functions_noblock',stringsAsFactors=F) %>%
  split(., seq(nrow(.)))
fields.cog.toprocess = pblapply(field.input,FUN=extract_pheno,dat.dic=fields.cog.outcome)

for (d in fields.cog.toprocess){
  if(colnames(d$phenotype)[1]==colnames(fields.cog.toprocess[[1]]$phenotype)[1]){
    dat.cog=d$phenotype[,c(2,1)]
    count.n=d$count
    count.n$type=rownames(count.n)
  }else{
    tmp.pheno=d$phenotype
    dat.cog=merge(dat.cog,tmp.pheno,by='f.eid',all.x=T)
    tmp.count.n=d$count
    tmp.count.n$type=rownames(tmp.count.n)
    count.n=merge(count.n,tmp.count.n,by='type',all.x=T)
  }
}

#saveRDS(count.n,file='data/countn.cognition.rds')

# calculate trail making outcomes ----------------------------------------------------------------------------

dat.cog$f.20157.diffs=dat.cog$f.20157-dat.cog$f.20156
dat.cog$f.20248.diffs=dat.cog$f.20248-dat.cog$f.20247
dat.cog$f.20155.diffs=dat.cog$f.20155-dat.cog$f.20149
dat.cog$f.20148.diffs=dat.cog$f.20148-dat.cog$f.20147

ls.trailmaking = paste0('f.',fields.cog.outcome$FieldID[fields.cog.outcome$Path=='Online follow-up > Cognitive function online > Trail making'])
dat.cog = dat.cog[,!(colnames(dat.cog) %in% ls.trailmaking)]

saveRDS(dat.cog,file=f.output_data)

# Data dictionary ----------------------------------------------------------------------------
fields.cog.new = fields.cog.outcome[fields.cog.outcome$Path!='Online follow-up > Cognitive function online > Trail making',]
fields.cog.new = rbind(fields.cog.new,fields.cog.outcome[grep('20157|20248|20155|20148',fields.cog.outcome$FieldID),])

fields.cog.new$category='Cognition'
fields.cog.new$field_tag=fields.cog.new$FieldID
fields.cog.new$field_tag[fields.cog.new$FieldID %in% c(20157,20248,20155,20148)]=paste0(fields.cog.new$field_tag[fields.cog.new$FieldID %in% c(20157,20248,20155,20148)],'.diffs')
fields.cog.new$field_used=fields.cog.new$FieldID
fields.cog.new$field_used[fields.cog.new$FieldID %in% c(20157,20248,20155,20148)]=c('20157,20156','20248,20247','20155,20149','20148,20147')

write.table(fields.cog.new,file=f.output_dictionary,sep='\t',quote=F,row.names=F,col.names=T)
