setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/')
library(dplyr)
library(pbapply)

# Load data    -----------------------------------------------------------------------------------------------
# Data dictionary
dat.dic=read.csv('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/Data_Dictionary_Showcase.csv',header=T,stringsAsFactors=F)
field.by.rls=readRDS('data/data_dictionary/fields_all.rds')
fields.all = merge(field.by.rls,dat.dic,by.x='field_name',by.y='FieldID',all.x=T)

fields.diet = fields.all[grep('Diet',fields.all$Path),]

write.table(fields.diet,file='data/data_dictionary/fields.diet_for_process.txt',row.names=F,col.names=T,sep='\t',quote=T)

# Data

diet=readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/2020-10-imaging-ukb40531/Touchscreen.rds')
diet.online=readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/2018-10-phenotypes-ukb24262/DietByHourRecall.rds')


# Process data fields    -------------------------------------------------------------------------------------

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

fields.diet = fields.diet %>%
                remove_category(.,kw='(pilot)',col.n='Field') %>%
                remove_category(.,kw='Reason',col.n='Field') %>%
                (function(x) x[,!colnames(x) %in% c(100020,20085,20086)]) %>%
                remove_category(.,kw='type',col.n='Field') %>%
                remove_category(.,kw='Type',col.n='Field')
                



                