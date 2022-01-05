library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.fields = args[1] 
d.pheno = args[2]
d.mdd = args[3]
f.output = args[4]
f.mismatch = args[5]

# f.fields = 'data/MR/MR_pheno_ls_UKB.RData'
# d.mdd = 'data/MR/MR_InterFiles_topSNP/mdd_exposure_outcome'
# d.pheno = 'data/MR/MR_InterFiles_topSNP/pheno_exposure_outcome'
# f.output = 'data/MR/ls.mr_analysis_topSNP.rds'
# f.mismatch = 'data/MR/ls.mr_mismatch_topSNP.RData'


# Load field info ---------------------------------------------------------

load(f.fields)

field.info = fields.all %>% select(field_tag,Field,category)

# Check matching phenotypes -----------------------------------------------

get_name <- function(x){
  x.out = x %>% basename %>% 
    gsub('.MDD.outcome_dat','',.,fixed = T) %>% 
    gsub('.outcome_dat','',.,fixed = T) %>% 
    gsub('.exposure_dat','',.,fixed = T)
  return(x.out)
}

ls.pheno_exposure = list.files(path = d.pheno, pattern='exposure_dat',full.names = T) %>% get_name
ls.pheno_outcome = list.files(path = d.pheno, pattern='outcome_dat',full.names = T) %>% get_name
ls.mdd_outcome = list.files(path = d.mdd, pattern='outcome_dat',full.names = T) %>% get_name
ls.mdd_exposure = list.files(path = d.mdd, pattern='exposure_dat',full.names = T) %>% get_name


ls.pheno_valid = Reduce(intersect, list(ls.pheno_exposure,ls.pheno_outcome,ls.mdd_outcome))

ls.mismatch.pheno_exposure = ls.pheno_exposure %>% .[!. %in% ls.pheno_valid]
ls.mismatch.pheno_outcome = ls.pheno_outcome %>% .[!. %in% ls.pheno_valid]
ls.mismatch.mdd_outcome = ls.mdd_outcome %>% .[!. %in% ls.pheno_valid]

save(ls.mismatch.mdd_outcome,ls.mismatch.pheno_exposure,
     ls.mismatch.pheno_outcome,file=f.mismatch)


# Create inputs for MR analysis -------------------------------------------

f.pheno_exposure = list.files(path = d.pheno, pattern='exposure_dat',full.names = T) 
f.pheno_outcome = list.files(path = d.pheno, pattern='outcome_dat',full.names = T) 
f.mdd_outcome = list.files(path = d.mdd, pattern='outcome_dat',full.names = T)
f.mdd_exposure = list.files(path = d.mdd, pattern='exposure_dat',full.names = T) 


inputs.mr.mdd_pheno = data.frame(file_exposure = f.mdd_exposure,
                                 file_outcome = f.pheno_outcome,
                                 stringsAsFactors = F) %>% 
  mutate(name_exposure=get_name(file_exposure),
         name_outcome=get_name(file_outcome)) %>% 
  mutate(label_exposure=name_exposure) %>% 
  left_join(field.info,by=c('name_outcome'='field_tag')) %>% 
  rename(label_outcome=Field) %>% 
  mutate(output_prefix=paste0(name_exposure,'_TO_',name_outcome)) %>% 
  select(file_exposure,file_outcome,
         name_exposure,name_outcome,
         label_exposure,label_outcome,
         category,output_prefix)

inputs.mr.pheno_mr = data.frame(file_exposure = f.pheno_exposure,
                                 file_outcome = f.mdd_outcome,
                                 stringsAsFactors = F) %>% 
  mutate(name_exposure=get_name(file_exposure),
         name_outcome='MDD') %>% 
  mutate(label_outcome=name_outcome) %>% 
  left_join(field.info,by=c('name_exposure'='field_tag')) %>% 
  rename(label_exposure=Field) %>% 
  mutate(output_prefix=paste0(name_exposure,'_TO_',name_outcome)) %>% 
  select(file_exposure,file_outcome,
         name_exposure,name_outcome,
         label_exposure,label_outcome,
         category,output_prefix)

inputs.all = rbind(inputs.mr.mdd_pheno,inputs.mr.pheno_mr)

saveRDS(inputs.all,file=f.output)