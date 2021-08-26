library(dplyr)
library(data.table)
library(readr)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dictionary = args[1]             # 
f.brain_category = args[2]         # 
f.IDP = args[3]                    # 
f.IDP_qccov = args[4]              # 
f.recruit = args[5]                # 
f.output_data = args[6]
f.output_dictionary = args[7]

# f.dictionary = 'data/Data_Dictionary_Showcase.csv'
# f.brain_category = 'results/phewas/data_dictionary/fields.imaging_phenotype.txt'
# f.IDP = 'data/2021-04-phenotypes-ukb44797/Imaging.rds'
# f.IDP_qccov = 'results/phewas/data_dictionary/fields.brain_imaging_QC_cov.txt'
# f.recruit = 'data/2021-04-phenotypes-ukb44797/Recruitment.rds'


# Load data ---------------------------------------------------------------

fields.all = fread(f.dictionary,header=T,stringsAsFactors=F) %>% data.frame
recruit = readRDS(f.recruit) %>% 
  select(f.eid,f.53.2.0,f.53.3.0)
# remove those who don't have a date for imaging assessments
IDP = readRDS(f.IDP) %>% 
  merge(.,recruit,by='f.eid') %>% 
  filter(!is.na(f.53.2.0)|!is.na(f.53.3.0)) 

# Select phenotypes to process --------------------------------------------

# Brain imaging data
remove_category <- function(x,kw,col.n,exact=F){
  if (exact==T){
    new.x=x[!grepl(paste0('^',kw,'$'),x[,col.n]),]
  }else{
    new.x=x[!grepl(kw,x[,col.n]),]
  }
  return(new.x)
}

fields.IDP <- fields.all[grep('Imaging > Brain MRI',fields.all$Path),] %>%
  remove_category(.,kw='Bulk',col.n='ItemType') %>% 
  remove_category(.,kw='L1 in tract|L2 in tract|L3 in tract|Mean L1|Mean L2|Mean L3',col.n='Field') %>% 
  remove_category(.,kw='Eprime|eprime',col.n='Field') %>% 
  .[!.$Coding %in% 1317,]
  

# Brain Imaging phenotypes (data dictionary and data)   ----------------------------------------------------------------

# Keep available fields in the data dictionary
ls.fields.available = colnames(IDP) %>% .[!. %in% 'f.eid'] %>% 
  strsplit(., '\\.') %>% 
  bind_cols %>% 
  t %>% data.frame(.,stringsAsFactors=F) %>% 
  .$X2 %>% as.numeric

fields.IDP = fields.IDP %>% 
  .[.$FieldID %in% ls.fields.available,]

# Only choose ASEG, DesikanPial, Desikan GW, DesikanFS, BAexvivo and bulk_tissue (others are duplicate measures using different atlases)
# for cortical and subcortical measures.
# Data: IDP.phenotype
# Dictionary: fields.category.process

ls.imaging.categories = read.delim(f.brain_category,header=T,stringsAsFactors=F)

add_category<-function(kw,col,rplc,cate.dat){
  cate.dat$category[grep(kw,cate.dat[,col])]=rplc
  return(cate.dat)
}
fields.category.process=fields.IDP %>% 
  .[.$Category %in% ls.imaging.categories$category,]
fields.category.process$category=NA
for (i in 1:nrow(ls.imaging.categories)){
  tmp.no.cate=ls.imaging.categories$category[i]
  tmp.label=ls.imaging.categories$label[i]
  
  fields.category.process <- add_category(tmp.no.cate,'Category',tmp.label,cate.dat=fields.category.process)
}

tmp.FieldID.tokeep=paste0('f.',fields.category.process$FieldID,'.2.0')
IDP.phenotype = IDP %>% 
  .[,colnames(.) %in% c('f.eid',tmp.FieldID.tokeep)]

# Combine left and right hemisphere
# Data: IDP.lnr
# Dictionary: fields.category.process.lnr

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

left.ref=fields.category.process[grep('left',fields.category.process$Field),]

right.ref = left.ref %>% 
  split(.,seq(nrow(.))) %>% 
  lapply(.,find_right,tmp.ref=fields.category.process) %>% 
  bind_rows

left.fields=left.ref$FieldID
right.fields=right.ref$FieldID

targetdata=IDP.phenotype
left.dat=targetdata[,paste0('f.',left.fields,'.2.0')]
right.dat=targetdata[,paste0('f.',right.fields,'.2.0')]

for (i in 1:ncol(left.dat)){
  aveg.source=cbind(left.dat[,i],right.dat[,i])
  tmp.aveg.dat=data.frame(f.eid=targetdata$f.eid,rowMeans(aveg.source,na.rm=T))
  tmp.aveg.dat[rowSums(!is.na(aveg.source))==0,2]=NA
  colnames(tmp.aveg.dat)[2]=gsub('\\.2\\.0','\\.lnr',colnames(left.dat)[i])
  if (i==1){
    aveg.dat=tmp.aveg.dat
  }else{
    aveg.dat=merge(aveg.dat,tmp.aveg.dat,by='f.eid',all.x=T)
  }
}

rm.fields=c(right.fields,left.fields)
rm.fields=paste0('f.',rm.fields,'.2.0')

IDP.lnr=IDP.phenotype[,!colnames(IDP.phenotype) %in% rm.fields] %>% 
    merge(.,aveg.dat,by='f.eid',all.x=T)

fields.category.process = fields.category.process %>% 
  .[!.$FieldID %in% right.fields,] %>% 
  mutate(field_tag=paste0('f.',FieldID),field_used=as.character(FieldID))

lnr.replacement.input = data.frame(FieldID=left.fields,
                                   used=paste0(left.fields,',',right.fields),
                                   stringsAsFactors=F) 

fields.category.process.lnr = fields.category.process
fields.category.process.lnr[match(left.fields,fields.category.process.lnr$FieldID),'field_used']=
  paste0(left.fields,',',right.fields)
fields.category.process.lnr[match(left.fields,fields.category.process.lnr$FieldID),'field_tag']=
  fields.category.process.lnr[match(left.fields,fields.category.process.lnr$FieldID),'field_tag'] %>%
  paste0(.,'.lnr')



# Brain Imaging QC and covariates (data dictionary and data)   ---------------------------------------------------------
# QC cov items
fields.QC_cov=read.delim(f.IDP_qccov,header=T,sep=';',stringsAsFactors=F)
fields.QC_cov=fields.all %>% 
  .[.$FieldID %in% fields.QC_cov$FieldID,] %>% 
  .[.$FieldID %in% ls.fields.available,] %>% 
  mutate(category='Brain imaging QC and covariates',
         field_tag=paste0('f.',FieldID,'.2.0'),
         field_used=as.character(FieldID))

IDP.cov=IDP %>% .[,colnames(.) %in% c('f.eid',paste0('f.',fields.QC_cov$FieldID,'.2.0'))]



# Save data and dictionary ------------------------------------------------

# Data
IDP.lnr=IDP.lnr %>% 
  .[,! colnames(.) %in% colnames(IDP.cov)[!grepl('f.eid',colnames(IDP.cov))]]
IDP.processed=merge(IDP.cov,IDP.lnr,by='f.eid',all.x=T)

saveRDS(IDP.processed,f.output_data)

# Dictionary
fields.BrainImaging = fields.category.process.lnr %>% 
  .[!fields.category.process.lnr$FieldID %in% fields.QC_cov$FieldID,] %>% 
  rbind(.,fields.QC_cov)

update_tag <- function(x,tmp.tag){
   real.tag = colnames(x) %>%
    .[grep(paste0('^',tmp.tag,'\\.|^',tmp.tag,'$'),.)]
   return(real.tag)
}

tag.update = fields.BrainImaging$field_tag %>%
   as.list %>%
   lapply(.,FUN=update_tag,x=IDP.processed) %>%
   unlist %>%
   as.character

fields.BrainImaging$field_tag = tag.update

write.table(fields.BrainImaging,file=f.output_dictionary,sep='\t',quote=F,row.names=F,col.names=T)
