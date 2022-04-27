#!/usr/bin/Rscript
library(data.table)
gwas<-fread('results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.11.gz')

gwas<-gwas[,c('SNP','A1','A2','FRQ_A_525197','FRQ_U_3362335','OR','SE','P','Nco','Nca'), with=F]

gwas$FREQ<-((gwas$FRQ_A_525197 * gwas$Nca) + (gwas$FRQ_U_3362335 * gwas$Nco)) / (gwas$Nca + gwas$Nco)
gwas$FRQ_A_525197<-NULL
gwas$FRQ_U_3362335<-NULL
gwas$BETA<-log(gwas$OR)
gwas$OR<-NULL
gwas$N<-gwas$Nca+gwas$Nco

gwas<-gwas[,c('SNP','A1','A2','FREQ','BETA','SE','P','N'), with=F]
names(gwas)<-c('SNP','A1','A2','freq','b','se','p','N')

fwrite(gwas, 'results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.11_COJO.txt', quote=F, na='NA', sep=' ')

