library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.PRS.new = args[1]             # f.PRS.new = 'data/mdd_prs_imputed.all_score'
f.PRS.old = args[2]             # f.PRS.old = 'data/mdd_prs_imputed_oldGWAS.all_score'
f.output_data = args[3]

# Reformat PRS ------------------------------------------------------------
PRS.new = fread(f.PRS.new,header=T) %>% 
  select(f.eid=IID,starts_with('Pt'))
colnames(PRS.new) = colnames(PRS.new) %>%
     gsub('Pt_','MDD3_Pt_',.) %>%
     gsub('e-','_',.)

PRS.old = fread(f.PRS.old,header=T) %>% 
  select(f.eid=IID,starts_with('Pt'))
colnames(PRS.old) = colnames(PRS.old) %>% 
     gsub('Pt_','Howard_Pt_',.) %>%
     gsub('e-','_',.)

PRS = left_join(PRS.new,PRS.old,by='f.eid')

# Save data ---------------------------------------------------------------

saveRDS(PRS,file=f.output_data)
