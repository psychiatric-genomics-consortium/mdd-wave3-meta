library(dplyr)
library(data.table)
library(readr)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dictionary = args[1]             # f.dictionary = 'data/Data_Dictionary_Showcase.csv'
f.baseline = args[2]               # f.baseline = 'data/2021-04-phenotypes-ukb44797/BaselineCharacteristics.rds'
f.recruit = args[3]                # f.recruit = 'data/2021-04-phenotypes-ukb44797/Recruitment.rds'
f.genomic = args[4]                # f.genomic = 'data/2021-04-phenotypes-ukb44797/derived/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.covars.rds'
f.online = args[5]                 # f.online = 'data/2021-04-phenotypes-ukb44797/MentalHealth.rds'
f.output_data = args[6]
f.output_dictionary = args[7]

# Load data ---------------------------------------------------------------

ls.covs = c(21022,34,52,21003,20400,53,54,31)
fields.all = fread(f.dictionary,header=T,stringsAsFactors=F) %>% data.frame %>% 
  .[.$FieldID %in% ls.covs,]
covs = readRDS(f.baseline) %>% 
  left_join(readRDS(f.recruit)) %>% 
  left_join(readRDS(f.online) %>% select(f.eid,matches('f.20400'))) %>% 
  select(f.eid,matches(paste0('f.',ls.covs,collapse = '|'))) %>% 
  select(-matches('.1.0$|.3.0$'))

genomic.covs = readRDS(f.genomic)
colnames(genomic.covs)[2:ncol(genomic.covs)] = 
  paste0('cov.',colnames(genomic.covs)[2:ncol(genomic.covs)])

# Process covs ------------------------------------------------------------
# Data
covs.output = covs %>% 
  mutate(cov.age_instance0 = f.21003.0.0,cov.age_instance2 = f.21003.2.0,
         cov.age_mhq = ((f.20400.0.0-f.53.0.0)/365+f.21003.0.0) %>% as.numeric,
         cov.assessment_centre_instance0 = f.54.0.0 %>% as.factor,
         cov.assessment_centre_instance2 = f.54.2.0 %>% as.factor,
         cov.sex = f.31.0.0 %>% as.factor) %>% 
  select(f.eid,starts_with('cov.')) %>% 
  right_join(genomic.covs %>% select(-num_range("cov.pc", 11:40)),by='f.eid')

# Dictionary
fields.output = fields.all %>% 
  .[.$FieldID %in% c(21003,54),] %>% 
  rbind(.,.) %>% 
  rbind(.,fields.all[fields.all$FieldID==31,]) %>%
  mutate(category='covariates',
         field_tag=c('cov.assessment_centre_instance0',
                     'cov.age_instance0',
                     'cov.assessment_centre_instance2',
                     'cov.age_instance2',
                     'cov.sex'),
         field_used=as.character(FieldID))

fields.output[grep('instance0',fields.output$field_tag),'Notes'] = 
  paste0(fields.output[grep('instance0',fields.output$field_tag),'Notes'],
         ' (instance 0, baseline assessment)')

fields.output[grep('instance2',fields.output$field_tag),'Notes'] = 
  paste0(fields.output[grep('instance2',fields.output$field_tag),'Notes'],
         ' (instance 2, imaging assessment)')

fields.cov.derived = fields.all[1:(ncol(covs.output)-1-nrow(fields.output)),] %>% 
  mutate(FieldID=colnames(covs.output) %>% .[!. %in% c('f.eid',fields.output$field_tag)]) %>%
  mutate(Link='https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=263',
         Notes='Phenotypes are derived covariates for genetic data') %>% 
  mutate_at(vars(-FieldID,-Link),function(x) x=NA) %>% 
  mutate(Participants=colSums(!is.na(covs.output[,
                                                 !colnames(covs.output) %in% c('f.eid',fields.output$field_tag)]))) %>% 
  mutate(category='covariates',field_tag=Field,field_used=NA) %>% 
  mutate(Field=gsub('cov.','Covariate: ',FieldID) %>% 
           gsub('age_','Age for ',.) %>% 
           gsub('\\.',' ',.) %>% 
           gsub('pc','genetic PC',.) %>% 
           gsub('mhq','Online mental health questionnaire',.),
         field_tag=FieldID)
fields.cov.derived[fields.cov.derived$field_tag=='cov.age_mhq','field_used']=
  '20400,53,21003'
rownames(fields.cov.derived)=NULL

fields.output=rbind(fields.output,fields.cov.derived)

update_tag <- function(x,tmp.tag){
  real.tag = colnames(x) %>%
    .[grep(paste0('^',tmp.tag,'\\.|^',tmp.tag,'$'),.)]
  return(real.tag)
}

tag.update = fields.output$field_tag %>%
  as.list %>%
  lapply(.,FUN=update_tag,x=covs.output) %>%
  unlist %>%
  as.character

fields.output$field_tag = tag.update



# Save data and dictionary ------------------------------------------------

saveRDS(covs.output,file=f.output_data)
write.table(fields.output,
            file=f.output_dictionary,sep='\t',quote=F,row.names=F,col.names=T)

