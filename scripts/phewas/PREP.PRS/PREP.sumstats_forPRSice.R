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
  filter(FRQ_A_235237>0.005,FRQ_A_235237<0.995,FRQ_U_2332871>0.005,FRQ_U_2332871<0.995) %>% 
  filter(Neff_half>5000) %>% 
  filter(Nca/Nco>0.005)

# write summstats
write_tsv(summstats,fname.out,col_names = T)