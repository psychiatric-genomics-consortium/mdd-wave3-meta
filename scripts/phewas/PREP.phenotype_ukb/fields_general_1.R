library(dplyr)
library(data.table)
library(readr)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dictionary = args[1]
f.experiment_procedure = args[2]
f.covered_by_sumscores = args[3]
f.output = args[4]



# Load data ---------------------------------------------------------------
dat.dic=fread(f.dictionary,header=T,stringsAsFactors=F) %>% as.data.frame


# Create loose fields to process ------------------------------------------------

remove_category <- function(x,kw,col.n,exact=F){
  if (exact==T){
    new.x=x[!grepl(paste0('^',kw,'$'),x[,col.n]),]
  }else{
    new.x=x[!grepl(kw,x[,col.n]),]
  }
  return(new.x)
}

fields.all = dat.dic
# Processed data chunks: imaging, MHQ, cognitive functions, health-related outcomes, 24-hr diet recall and touchscreen diet, alcohol consumption
# Other non-phenotype variables: genomic matrices, notes, device ID, sampling notes and non-data categories
ls.categories.rm = c(2,152,100,103,105:112,123,127,129,134,135,148,150,
                     190:197,220:222,501:506,1009:1013,1101,1102,9081,18518,100002,
                     100004,100016,100023:100032,100052,100060:100061,100075,100081,
                     100085:100087,100095,100096) %>% 
  paste0('^',.,'$',collapse = '|')
fields.loose = fields.all %>%
  filter(.,!is.na(Field)) %>%
  remove_category(.,kw='Online follow-up > Mental health',col.n='Path') %>%
  remove_category(.,kw='UK Biobank Assessment Centre > Imaging > Brain MRI',col.n='Path') %>%
  remove_category(.,kw='UK Biobank Assessment Centre > Cognitive function',col.n='Path') %>%
  remove_category(.,kw='Online follow-up > Cognitive function',col.n='Path') %>%
  remove_category(.,kw='Health-related outcomes ',col.n='Path') %>%
  remove_category(.,kw='Online follow-up > Diet by 24-hour recall |Touchscreen > Lifestyle and environment > Alcohol',col.n='Path') %>%
  remove_category(.,kw='UK Biobank Assessment Centre > Touchscreen > Lifestyle and environment > Diet',col.n='Path') %>%
  remove_category(.,kw='UK Biobank Assessment Centre > Touchscreen > Lifestyle and environment > Alcohol',col.n='Path') %>%
  remove_category(.,kw='Online follow-up > Diet',col.n='Path') %>%
  remove_category(.,kw='Genomics',col.n='Path') %>%
  remove_category(.,kw='UK Biobank Assessment Centre > Procedural metrics',col.n='Path') %>%
  remove_category(.,kw='Reason for ',col.n='Field') %>%
  remove_category(.,kw='device ID|Identifier|error indicator|household related to participant',col.n='Field') %>%
  remove_category(.,kw='UK Biobank Assessment Centre > Biological sampling',col.n='Path') %>%
  remove_category(.,kw='(pilot)|Invitation| method|Which eye(s)|Method of ',col.n='Field') %>%
  remove_category(.,kw=' code|Duration of questionnaire|questionnaire completed',col.n='Field') %>%
  remove_category(.,kw='Bulk',col.n='ItemType') %>% 
  remove_category(.,kw='Which eye(s)',col.n='Path',exact = T) %>% 
  remove_category(.,kw=ls.categories.rm,col.n='Category')
  

rm.others=read.csv(f.experiment_procedure,header=T,sep='\t',stringsAsFactors=F) %>%
  rbind(.,read.csv(f.covered_by_sumscores,header=T,sep='\t',stringsAsFactors=F))

fields.loose=fields.loose[!fields.loose$FieldID %in% as.numeric(rm.others$field_name),]

write.table(fields.loose,file=f.output,quote=F,sep='\t',row.names = F,col.names = T)

