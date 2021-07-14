library(dplyr)
library(data.table)
library(readr)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

fam.file = args[1]
f.valid_subj = args[2]
output.file = args[3]

# Make dummy phenotype file -----------------------------------------------

fam.local = fread(fam.file,header = F,stringsAsFactors = F) 
valid.subj = fread(f.valid_subj,header=F,stringsAsFactors = F) %>% .$V1

dummy.pheno = fam.local %>% 
  select(FID = V1, IID = V2) %>% 
  .[.$IID %in% valid.subj,] %>% 
  mutate(dummy_pheno = runif(nrow(.),0,1)) %>% 
  mutate(dummy_pheno = round(dummy_pheno))

write_tsv(dummy.pheno,output.file,col_names = T)
