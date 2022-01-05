# XShen 13/07/2021
# Basic settings -----------------------------------------------------------------
library(dplyr)
library(data.table)
library(readr)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

option_list <- list(
  make_option('--sumstats', type='character', help="sumstats file", action='store'),
  make_option('--out', type='character', help="output file", action='store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

fname.sumstats = opt$sumstats
fname.out = opt$out

# Summary statistics  ---------------------------------------------------------------------
summstats=read.delim(fname.sumstats,stringsAsFactors=F)

# Reformat and QC   
summstats=filter(summstats,INFO>=0.9,!is.na(INFO)) %>% 
  filter(FRQ_A_470188>0.005,FRQ_A_470188<0.995,FRQ_U_2752545>0.005,FRQ_U_2752545<0.995) %>% 
  filter(Neff_half>(Neff_half*0.8)) %>% 
  filter(Nca/Nco>0.005)

# write summstats
write_tsv(summstats,fname.out,col_names = T)
summstats %>% 
  as.data.frame %>% 
  .[,'SNP'] %>% 
write.table(.,'data/mdd3.snps',col.names=F,row.names=F,quote=F,sep='\n')