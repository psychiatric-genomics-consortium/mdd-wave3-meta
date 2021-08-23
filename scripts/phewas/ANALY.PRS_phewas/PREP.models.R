library(dplyr)
library(data.table)
library(readr)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dat.brain = args[1]
f.dat.cog = args[2]
f.dat.diet = args[3]
f.dat.activity = args[4]
f.dat.mental = args[5]
f.dat.loose = args[6]
f.dat.cov = args[7]

f.dic_dir = args[8]             
f.PRS = args[9]                 
f.output_data = args[10]

# f.dat.brain = 'data/dat.imaging_chunk.rds'
# f.dat.cog = 'data/dat.cognition_chunk.rds'
# f.dat.diet = 'data/dat.diet_chunk.rds'
# f.dat.activity = 'data/dat.activity_chunk.rds'
# f.dat.mental = 'data/dat.mental_health_chunk.rds'
# f.dat.loose = 'data/dat.loose_field_chunk.rds'
# f.dat.cov = 'data/dat.addional_covariates_chunk.rds'
# f.dic_dir = 'results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt'
# f.PRS = 'data/PRS_all.rds'


# Load data ---------------------------------------------------------------

# Data
fs = ls(pattern = 'f.dat') %>% as.list %>% lapply(get)

phewas.dat = fs %>% lapply(readRDS) %>% 
  Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="f.eid",all=T), .) %>% 
  right_join(.,readRDS(f.PRS),by='f.eid')

#n.rm = colSums(!is.na(phewas.dat))<2000 
#ls.rm = colnames(phewas.dat)[n.rm] %>% .[!. %in% 'f.eid']
#phewas.dat = phewas.dat %>% 
#  .[,!colnames(.) %in% ls.rm]

PRS = readRDS(f.PRS)

# Data dictionary
dic.path = f.dic_dir %>% strsplit(.,'/') %>% unlist %>% 
  .[1:(length(.)-1)] %>% 
  paste0(.,collapse = '/')

fields.all = list.files(dic.path,pattern='^fields.final',full.names = T,include.dirs = T) %>% 
  lapply(fread,stringsAsFactors=F,header=T,sep='\t',quote="\"") %>%
  lapply(as.data.frame) %>% 
  lapply(mutate,field_used=as.character(field_used),FieldID=as.character(FieldID)) %>% 
  Reduce(function(dtf1,dtf2) bind_rows(dtf1,dtf2[,colnames(dtf1)]), .) #%>% 
  #.[!.$field_tag %in% ls.rm,]

fields.qc_cov = fields.all %>% 
  .[grepl('cov|Cov',.$category),]

fields.toi = fields.all %>% 
  .[!.$FieldID %in% fields.qc_cov$FieldID,] %>% 
  mutate(Field=gsub(' (left)','',Field)) %>% 
  mutate(Field=gsub('(left)','',Field))

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
saveRDS(fields.toi,file = 'data/fields_toi.rds')
