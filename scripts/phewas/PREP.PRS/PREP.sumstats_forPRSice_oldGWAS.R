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
ref=fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ucsc_annotation/hg19/snp151Common_chr_bp_rs.txt') %>% 
  rename(CHR=V1,BP=V2,SNP=V3) %>% 
  filter(CHR!='X',CHR!='Y') %>% 
  mutate(CHR=gsub('chr','',CHR) %>% as.numeric)

# Reformat and QC   
summstats=summstats %>% 
  filter(Freq1>0.005,Freq1<0.995) %>% 
  select(SNP=MarkerName,A1=Allele1,A2=Allele2,SE=StdErr,P=P.value,BETA=Effect) %>% 
  left_join(.,ref,by='SNP')

# write summstats
write_tsv(summstats,fname.out,col_names = T)