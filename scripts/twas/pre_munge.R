#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
                   
library(data.table)
ss<-fread('results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.11.gz')

ss<-ss[,c(1:5,8:11,19),with=F]
ss$N<-ss$Neff_half*2

fwrite(ss, 'results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.11_premunged.gz', sep=' ', na='NA', quote=F)
