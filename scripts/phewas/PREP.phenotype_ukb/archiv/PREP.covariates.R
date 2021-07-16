setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/')
library(dplyr)
library(pbapply)


recruit=readRDS('/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/2020-10-imaging-ukb40531/Recruitment.rds')
recruit=recruit[,c('f.eid','f.54.0.0','f.21003.0.0')]
colnames(recruit)[2:3]=c('assessment_centre','age')

baseline=readRDS('/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/2018-10-phenotypes-ukb24262/BaselineCharacteristics.rds')
baseline=baseline[,c('f.eid','f.31.0.0')]
colnames(baseline)[2]='sex'

covariates=merge(recruit,baseline,by='f.eid',all.x=T)

genetic.covs=readRDS('/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/Share/ukb_array_pcs_v2.rds')
genetic.covs=genetic.covs[,c(1:2,3:18)]

covariates=merge(covariates, genetic.covs,by='f.eid',all.x=T)


imaging.covs=readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/2020-10-imaging-ukb40531/Imaging.rds')
cov.fields = c(25731,25757:25759)
imaging.covs=imaging.covs[,c('f.eid',paste0('f.',cov.fields,'.2.0'))]
colnames(imaging.covs)[2:ncol(imaging.covs)]=c('est.ICV','pos.x','pos.y','pos.z')

covariates=merge(covariates,imaging.covs,by='f.eid',all.x=T)

ls.ID.include=read.delim('/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/Share/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id',header=F,sep=' ')
colnames(ls.ID.include)=c('f.eid','ID.tokeep')

covariates = covariates[covariates$f.eid %in% ls.ID.include$f.eid,]

saveRDS(covariates,file='data/covariates.rds')