library(dplyr)
library(data.table)
library(readr)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dic_dir = args[1]             # f.dic_dir = 'results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt'
f.PRS = args[2]                 # f.PRS = 'data/PRS_all.rds'
f.output_data = args[3]

# Load data ---------------------------------------------------------------

g.path = f.dic_dir %>% strsplit(.,'/') %>% unlist %>% 
  .[1:(length(.)-1)] %>% 
  paste0(.,collapse = '/')

fields.all = list.files(g.path,pattern='^fields.final',full.names = T,include.dirs = T) %>% 
  lapply(fread,stringsAsFactors=F,header=T,sep='\t',quote="\"") %>%
  lapply(as.data.frame) %>% 
  lapply(mutate,field_used=as.character(field_used),FieldID=as.character(FieldID)) %>% 
  Reduce(function(dtf1,dtf2) bind_rows(dtf1,dtf2[,colnames(dtf1)]), .)

fields.qc_cov = fields.all %>% 
  .[grepl('cov|Cov',.$category),]

fields.toi = fields.all %>% 
  .[!.$FieldID %in% fields.qc_cov$FieldID,] %>% 
  mutate(Field=gsub(' (left)','',Field)) %>% 
  mutate(Field=gsub('(left)','',Field))

PRS = readRDS(f.PRS)


# Build models ------------------------------------------------------------

##### dependent variables
ls.dep = fields.toi$field_tag

##### factors
ls.factor = colnames(PRS) %>% .[!.%in%'f.eid']

##### combine the two
ls.models=expand.grid(ls.dep,ls.factor,stringsAsFactors = F) %>% 
  as.data.frame(stringsAsFactors=F) %>% 
  rename(dependent=Var1,factor=Var2) %>% 
  mutate(covs='',est='glm') %>% 
  left_join(fields.toi[,c('field_tag','category')],by=c('dependent'='field_tag'))

##### covs
# Common covs
ls.models$covs=paste0('cov.pc',1:10) %>% 
  c(.,'cov.genotyping.array','cov.assessment_centre_instance0','cov.sex','cov.age_instance0') %>% 
  paste0(.,collapse = '+')

# Imaging phenotypes
ls.models$covs[grepl('atlas|brain|White matter',ls.models$dependent)]=
  ls.models$covs[grepl('atlas|brain|White matter',ls.models$dependent)] %>% 
  gsub('instance0','instance2',.) %>% 
  paste0(.,'+f.25756.2.0+f.25757.2.0+f.25758.2.0')

ls.models$covs[grepl('atlas',ls.models$dependent)]=
  ls.models$covs[grepl('atlas',ls.models$dependent)] %>% 
  paste0(.,'+f.25000.2.0') 

ls.models$covs[grepl('White matter',ls.models$dependent)]=
  ls.models$covs[grepl('White matter',ls.models$dependent)] %>% 
  paste0(.,'+f.25737.2.0') 

# Online questionnaire
ls.models$covs[grepl('Mental health',ls.models$dependent)&!grepl('^f.',ls.models$field_tag)]=
  ls.models$covs[grepl('Mental health',ls.models$dependent)&!grepl('^f.',ls.models$field_tag)] %>% 
  gsub('+cov.age_instance0','+cov.age_mhq',.)

# Physical activity
additional.covs.activity = fields.qc_cov$field_tag[grepl('Calibration coefficients|Wear duration',
                                                         fields.qc_cov$Field)] %>% 
  paste0(.,collapse = '+')
ls.models$covs[ls.models$dependent %in% c(1009,1010)]=
  ls.models$covs[ls.models$dependent %in% c(1009,1010)] %>% 
  paste0(.,'+',additional.covs.activity)


# Save model data ---------------------------------------------------------

saveRDS(ls.models,file = f.output_data)
