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
f.activity = args[3]                # f.activity = 'data/2021-04-phenotypes-ukb44797/PhysicalActivityMeasurement.rds'
f.output_data = args[4]
f.output_dictionary = args[5]


# Load data ---------------------------------------------------------------

fields.all=read.csv(f.dic,header=T,stringsAsFactors=F)
ref.coding=read.csv(f.coding,header=T,stringsAsFactors=F) 

# Derived data only, categories 1009 and 1010 included as data
# wearing duration (FieldID==90051) and calibration covariates (FieldID==90161:90169) included as covs
# Data quality after calibration (FieldID==90016) as QC variable ('Yes' included)

fields.activity.dat = fields.all %>% 
  .[grep('Additional exposures > Physical activity measurement',.$Path),] %>% 
  .[grep('Data',.$ItemType),] %>% 
  .[grep('Derived',.$Strata),] %>% 
  .[grep('1009|1010',.$Category),]
fields.activity.covs = fields.all %>% 
  .[grep(paste0(c(90051,90161:90169),collapse = '|'),.$FieldID),] 
fields.activity = rbind(fields.activity.dat,fields.activity.covs)

# Only include participants with good data quality after calibration
activity=readRDS(f.activity) %>% 
  filter(f.90016.0.0=='Yes') 

ls.fields.available = colnames(activity) %>% .[!. %in% 'f.eid'] %>% 
  strsplit(., '\\.') %>% 
  bind_cols %>% 
  t %>% data.frame(.,stringsAsFactors=F) %>% 
  .$X2 %>% as.numeric

fields.activity = fields.activity %>% 
  .[.$FieldID %in% ls.fields.available,]

# Note: all columns are numerical

# Process data and covariates ---------------------------------------------

activity.output = activity %>% 
  .[,c('f.eid',paste0('f.',fields.activity$FieldID,'.0.0'))]

# Process data dictionary -------------------------------------------------

field.activity.output = fields.activity
field.activity.output$category = 'Physical activity'
field.activity.output$category[field.activity.output$FieldID %in% fields.activity.covs$FieldID] = 
  'Physical activity QC and covariates'
field.activity.output = field.activity.output %>% 
  mutate(field_tag = paste0('f.',FieldID,'.0.0'), field_used = as.character(FieldID))

update_tag <- function(x,tmp.tag){
  real.tag = colnames(x) %>%
    .[grep(paste0('^',tmp.tag,'\\.|^',tmp.tag,'$'),.)]
  return(real.tag)
}

tag.update = field.activity.output$field_tag %>%
  as.list %>%
  lapply(.,FUN=update_tag,x=activity.output) %>%
  unlist %>%
  as.character

field.activity.output$field_tag = tag.update


# Save data, covariates and data dictionary -------------------------------

saveRDS(activity.output,file=f.output_data)
write.table(field.activity.output,
            file=f.output_dictionary,sep='\t',quote=T,row.names=F,col.names=T)

