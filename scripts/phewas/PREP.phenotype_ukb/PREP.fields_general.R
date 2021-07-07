setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/PGCMDD3-meta-PheWAS')

library(dplyr)

dat.dic=read.csv('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/Data_Dictionary_Showcase.csv',header=T,stringsAsFactors=F)

# Load fields by release   -----------------------------------------------------------------------------------

g.path = '/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/'

dirs = list.dirs(path=g.path,recursive=F)%>% (function(x) x[grep('-ukb',x)])
rls <- dirs %>%
        gsub(g.path,'',.) %>%
        strsplit(.,'-') %>%
        lapply(function(x) x[length(x)]) %>%
        rev(.) # rev to start from the newest release

get_field <- function(rls.input,dirs.input){
    tmp.dir = dirs.input[grep(rls.input,dirs.input)]
    field.file = paste0(tmp.dir,'/fields.ukb')
    fields = read.table(field.file,header=F)
    fields.rls = data.frame(rls=rls.input,field_name=fields$V1)
    return(fields.rls)
}

field.by.rls = lapply(X=rls,FUN=get_field,dirs.input=dirs) %>%
                bind_rows(.) %>%
                (function(x) x[!duplicated(x$field_name),])
saveRDS(field.by.rls,file='data/data_dictionary/fields_all.rds')   


# Remove data that need separate processing and the non-phenotype variables ---------------------------------------

remove_category <- function(x,kw,col.n,exact=F){
      if (exact==T){
          new.x=x[!grepl(paste0('^',kw,'$'),x[,col.n]),]
      }else{
      new.x=x[!grepl(kw,x[,col.n]),]
      }
      return(new.x)
}

fields.all = merge(field.by.rls,dat.dic,by.x='field_name',by.y='FieldID',all.x=T)
# Processed data chunks: imaging, MHQ, cognitive functions, health-related outcomes, 24-hr diet recall and touchscreen diet, alcohol consumption
# Other non-phenotype variables: genomic matrices, notes, device ID, sampling notes
fields.loose = fields.all %>%
                filter(.,!is.na(Field)) %>%
                remove_category(.,kw='Online follow-up > Mental health',col.n='Path') %>%
                remove_category(.,kw='UK Biobank Assessment Centre > Imaging > Brain MRI',col.n='Path') %>%
                remove_category(.,kw='UK Biobank Assessment Centre > Cognitive function',col.n='Path') %>%
                remove_category(.,kw='Online follow-up > Cognitive function',col.n='Path') %>%
                remove_category(.,kw='Health-related outcomes ',col.n='Path') %>%
                remove_category(.,kw='Online follow-up > Diet by 24-hour recall ',col.n='Path') %>%
                remove_category(.,kw='UK Biobank Assessment Centre > Touchscreen > Lifestyle and environment > Diet',col.n='Path') %>%
                remove_category(.,kw='UK Biobank Assessment Centre > Touchscreen > Lifestyle and environment > Alcohol',col.n='Path') %>%
                remove_category(.,kw='Genomics',col.n='Path') %>%
                remove_category(.,kw='UK Biobank Assessment Centre > Procedural metrics',col.n='Path') %>%
                remove_category(.,kw='Reason for ',col.n='Field') %>%
                remove_category(.,kw='device ID',col.n='Field') %>%
                remove_category(.,kw='UK Biobank Assessment Centre > Biological sampling',col.n='Path') %>%
                remove_category(.,kw='(pilot)|Invitation| method|Which eye(s)|Method of ',col.n='Field') %>%
                remove_category(.,kw=' code',col.n='Field') %>%
                remove_category(.,kw='Bulk',col.n='ItemType')                
                
write.table(fields.loose,file='data/data_dictionary/fields_for_manual_selection.csv',quote=F,col.names=T,row.names=F,sep='\t')

# The following two files were generated manually
# field_experiment_procedure.txt: text notes or variables for experiment procedures
# field_covered_by_summary_scores.txt: item scores that can be covered by general scores

rm.others=read.csv('data/data_dictionary/field_experiment_procedure.txt',header=T,sep='\t',stringsAsFactors=F)
rm.others=rbind(rm.others,read.csv('data/data_dictionary/field_covered_by_summary_scores.txt',header=T,sep='\t',stringsAsFactors=F))

fields.loose=fields.loose[!fields.loose$field_name %in% rm.others$field_name,]

write.table(fields.loose,file='data/data_dictionary/fields_loose_for_process.csv',quote=F,col.names=T,row.names=F,sep='\t')

# Run PREP.fields_separate_process.R after running this script