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

# f.cog_ac = 'data/2021-04-phenotypes-ukb44797/CognitiveFunction.rds'
# f.cog_online = 'data/2021-04-phenotypes-ukb44797/CognitiveFunctionOnline.rds'
# f.cogvar_QCed = 'results/phewas/data_dictionary/fields.cog_after_manual_check.txt'


# Load data ---------------------------------------------------------------

cognition=readRDS(f.cog_ac) %>% 
  left_join(.,readRDS(f.cog_online),by='f.eid')

ls.fields.available = colnames(cognition) %>% .[!. %in% 'f.eid'] %>% 
  strsplit(., '\\.') %>% 
  bind_cols %>% 
  t %>% data.frame(.,stringsAsFactors=F) %>% 
  .$X2 %>% as.numeric

fields.cog=read.delim(f.cogvar_QCed,header=T,stringsAsFactors=F) %>% 
  .[.$FieldID %in% ls.fields.available,] %>%
  select(-rls)


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
  select(matches(paste(paste0('f.',fields.cog$FieldID), collapse="|")))

# find unique list of task * instances
ls.task_instance = colnames(cog_functions_raw) %>% 
  .[ !. %in% 'f.eid'] %>% 
  strsplit(., '\\.') %>% 
  bind_cols %>% 
  t
rownames(ls.task_instance)=NULL
ls.fields_findblock <- apply(ls.task_instance[,1:3],1,paste,collapse=".") %>% 
  unique %>% paste0(.,'.')

# collapse block data: extract mean
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
field.input=data.frame(field.id=fields.cog$FieldID,dat='cog_functions_noblock',stringsAsFactors=F) %>%
  split(., seq(nrow(.)))
dat.cog.singleinstance = pblapply(field.input,FUN=extract_pheno,dat.dic=fields.cog) %>% 
  pblapply(.,function(x) x$phenotype) %>% 
  Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="f.eid"), .)

# calculate trail making outcomes ----------------------------------------------------------------------------

# Data
dat.cog.singleinstance = dat.cog.singleinstance %>% 
  mutate(f.20157.diffs = f.20157-f.20156,
         f.20248.diffs = f.20248-f.20247,
         f.20155.diffs = f.20155-f.20149,
         f.20148.diffs = f.20148-f.20147)

ls.trailmaking = paste0('f.',fields.cog$FieldID[fields.cog$Path=='Online follow-up > Cognitive function online > Trail making'])
dat.cog.singleinstance = dat.cog.singleinstance %>% 
  .[,!(colnames(.) %in% ls.trailmaking)]

# Update data dictionary
fields.cog.newTrailMaking = fields.cog %>% 
  .[.$Path!='Online follow-up > Cognitive function online > Trail making',] %>% 
  rbind(.,fields.cog[grep('20157|20248|20155|20148',fields.cog$FieldID),])

fields.cog.newTrailMaking =fields.cog.newTrailMaking %>% 
  mutate(category='Cognition',field_tag=paste0('f.',FieldID),field_used=FieldID)
fields.cog.newTrailMaking$field_tag[fields.cog.newTrailMaking$FieldID %in% c(20157,20248,20155,20148)]=
  paste0(fields.cog.newTrailMaking$field_tag[fields.cog.newTrailMaking$FieldID %in% c(20157,20248,20155,20148)],'.diffs')
fields.cog.newTrailMaking$field_used[fields.cog.newTrailMaking$FieldID %in% c(20157,20248,20155,20148)]=
  c('20157,20156','20248,20247','20155,20149','20148,20147')

fields.cog.newTrailMaking$Field[fields.cog.newTrailMaking$FieldID %in% c(20157,20248,20155,20148)]=
  gsub('trail #2','trail2-trail1',
       fields.cog.newTrailMaking$Field[fields.cog.newTrailMaking$FieldID %in% c(20157,20248,20155,20148)])

fields.cog.newTrailMaking = fields.cog.newTrailMaking %>% 
  mutate(field_used=as.character(field_used))
# Update field tag
update_tag <- function(x,tmp.tag){
  real.tag = colnames(x) %>%
    .[grep(paste0('^',tmp.tag,'\\.|^',tmp.tag,'$'),.)]
  return(real.tag)
}

tag.update = fields.cog.newTrailMaking$field_tag %>%
  as.list %>%
  lapply(.,FUN=update_tag,x=dat.cog.singleinstance) %>%
  unlist %>%
  as.character

fields.cog.newTrailMaking$field_tag = tag.update


# Save data and dictionary ------------------------------------------------

saveRDS(dat.cog.singleinstance,file=f.output_data)
write.table(fields.cog.newTrailMaking,file=f.output_dictionary,sep='\t',quote=T,row.names=F,col.names=T)
