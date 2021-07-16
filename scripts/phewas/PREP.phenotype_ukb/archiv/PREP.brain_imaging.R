setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/')
library(dplyr)
library(pbapply)

# Load data    -----------------------------------------------------------------------------------------------
# Data dictionary
dat.dic=read.csv('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/Data_Dictionary_Showcase.csv',header=T,stringsAsFactors=F)
field.by.rls=readRDS('data/data_dictionary/fields_all.rds')
fields.all = merge(field.by.rls,dat.dic,by.x='field_name',by.y='FieldID',all.x=T)

# Data
IDP=readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/2020-10-imaging-ukb40531/Imaging.rds')
recruit=readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/2020-10-imaging-ukb40531/Recruitment.rds')
IDP=filter(IDP,!is.na(recruit$f.53.2.0)|!is.na(recruit$f.53.3.0)) # remove people who don't have date of attendance at imaging assessments
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
              remove_category(.,kw='Bulk',col.n='ItemType') 
              
              
              
# Brain Imaging phenotypes (data dictionary and data)   ----------------------------------------------------------------

# Only choose ASEG, DesikanPial, Desikan GW, DesikanFS, BAexvivo and bulk_tissue (others are duplicate measures using different atlases)
# This applies only to cortical and subcortical measures, as multiple atlases were used.

ls.imaging.categories = read.delim('data/data_dictionary/fields.imaging_phenotype.txt',header=T,stringsAsFactors=F)

add_category<-function(kw,col,rplc,cate.dat){
    cate.dat$category[grep(kw,cate.dat[,col])]=rplc
    return(cate.dat)
}
fields.category.process=fields.IDP
fields.category.process$category=NA
for (i in 1:nrow(ls.imaging.categories)){
    tmp.no.cate=ls.imaging.categories$category[i]
    tmp.label=ls.imaging.categories$label[i]
    
    fields.category.process <- add_category(tmp.no.cate,'Category',tmp.label,cate.dat=fields.category.process)
}
fields.category.process=filter(fields.category.process,!is.na(category))

tmp.field_name.tokeep=c(paste0('f.',fields.category.process$field_name,'.2.0'),paste0('f.',fields.category.process$field_name,'.2.0'))
IDP.phenotype = IDP[,colnames(IDP) %in% c('f.eid',tmp.field_name.tokeep)]

# Combine left and right hemisphere

fields.category.process$Field=as.character(fields.category.process$Field)
left.fields=fields.category.process$field_name[grep('left',fields.category.process$Field)]
right.fields<- gsub('left','right',x=fields.category.process$Field[grep('left',fields.category.process$Field)]) %>%
                match(.,fields.category.process$Field) %>%
                (function(x) fields.category.process$field_name[x])
targetdata=IDP.phenotype
left.dat=targetdata[,paste0('f.',left.fields,'.2.0')]
right.dat=targetdata[,paste0('f.',right.fields,'.2.0')]


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
rm.fields=paste0('f.',rm.fields,'.2.0')

IDP.lnr=IDP.phenotype[,!colnames(IDP.phenotype) %in% rm.fields]
IDP.lnr=merge(IDP.lnr,aveg.dat,by='f.eid',all.x=T)

# Brain Imaging QC and covariates (data dictionary and data)   ---------------------------------------------------------
# QC cov items
fields.QC_cov=read.delim('data/data_dictionary/fields.brain_imaging_QC_cov.txt',header=T,sep=';',stringsAsFactors=F)
fields.QC_cov=fields.all[fields.all$field_name %in% fields.QC_cov$field_name,]
fields.QC_cov$category='Brain imaging QC and covariates'

IDP.cov=IDP[,colnames(IDP) %in% paste0('f.',fields.QC_cov$field_name,'.2.0')]
IDP.cov$f.eid=IDP$f.eid

# Save the new IDP   ---------------------------------------------------------
IDP.lnr=IDP.lnr[,!colnames(IDP.lnr) %in% colnames(IDP.cov)[2:ncol(IDP.cov)]]
IDP.processed=merge(IDP.cov,IDP.lnr,by='f.eid',all.x=T)

saveRDS(IDP.processed,'data/dat.imaging_chunk.rds')

# Update data dictionary   ---------------------------------------------------------
tmp.fields.block.process=rbind(fields.QC_cov,fields.category.process)
tmp.fields.block.process=tmp.fields.block.process[!duplicated(tmp.fields.block.process$field_name),]

rm.fields.dic <- right.fields
fields.imaging_lnr = tmp.fields.block.process[! tmp.fields.block.process$field_name %in% rm.fields.dic,]
fields.imaging_lnr$field_tag = as.character(fields.imaging_lnr$field_name)
fields.imaging_lnr$field_tag[fields.imaging_lnr$field_name %in% left.fields]=
      paste0(fields.imaging_lnr$field_tag[fields.imaging_lnr$field_name %in% left.fields],'2.0.lnr')
fields.imaging_lnr$fields_used=fields.imaging_lnr$field_tag
fields.imaging_lnr$fields_used[fields.imaging_lnr$field_name %in% left.fields]=
      paste0(fields.imaging_lnr$fields_used[fields.imaging_lnr$field_name %in% left.fields],
              fields.imaging_lnr$fields_used[fields.imaging_lnr$field_name %in% right.fields])

write.table(fields.imaging_lnr,file='data/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt',sep='\t',row.names=F,col.names=T)

