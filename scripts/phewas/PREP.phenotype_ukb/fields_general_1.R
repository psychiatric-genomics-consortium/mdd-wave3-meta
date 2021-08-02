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
dat.dic=read.csv(f.dictionary,header=T,stringsAsFactors=F)


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
  remove_category(.,kw='Online follow-up > Diet',col.n='Path') %>%
  remove_category(.,kw='Genomics',col.n='Path') %>%
  remove_category(.,kw='UK Biobank Assessment Centre > Procedural metrics',col.n='Path') %>%
  remove_category(.,kw='Reason for ',col.n='Field') %>%
  remove_category(.,kw='device ID',col.n='Field') %>%
  remove_category(.,kw='UK Biobank Assessment Centre > Biological sampling',col.n='Path') %>%
  remove_category(.,kw='(pilot)|Invitation| method|Which eye(s)|Method of ',col.n='Field') %>%
  remove_category(.,kw=' code',col.n='Field') %>%
  remove_category(.,kw='Bulk',col.n='ItemType')                

rm.others=read.csv(f.experiment_procedure,header=T,sep='\t',stringsAsFactors=F) %>% 
  rbind(.,read.csv(f.covered_by_sumscores,header=T,sep='\t',stringsAsFactors=F))

fields.loose=fields.loose[!fields.loose$FieldID %in% rm.others$FieldID,]

write.table(fields.loose,file=f.output,quote=F,col.names=T,row.names=F,sep='\t')

