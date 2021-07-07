# XShen 22/08/2020
# Basic settings -----------------------------------------------------------------
library(dplyr)

setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/')

# Summary statistics  ---------------------------------------------------------------------
summstats=read.delim('data/daner_pgc_mdd_noUKBB_eur_hg19_v3.229.09.txt',stringsAsFactors=F)

# Reformat and QC   
summstats=filter(summstats,INFO>=0.9,!is.na(INFO)) # Info > 0.9
# Max Nca+Nco=9302785 ?? possible?
summstats=filter(summstats,FRQ_A_176397>0.005,FRQ_A_176397<0.995) # MAF<0.005 
summstats=filter(summstats,FRQ_U_2064500>0.005,FRQ_U_2064500<0.995) # MAF<0.005

write.table(summstats,file='data/daner_pgc_mdd_noUKBB_eur_hg19_v3.229.09_forPRS.txt',sep=' ',quote=F,row.names=F,col.names=T)

# MDD phenotype for unrelated subj in UKBB ------------------------------------------------
mdd_touchscreen=readRDS('/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/mdd_pipeline_MAdams/ukb8238-mdd.mdd_phenotypes.rds')
mdd_CIDI<-readRDS('/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/mhq/MHQ.1907.ukb24262.Derived.rds') %>% (function(x) x[,c('f.eid','Depressed.Ever')])
        
mdd.phenotypes=merge(mdd_touchscreen,mdd_CIDI,by='f.eid',all.x=T)
mdd.phenotypes$IID=mdd.phenotypes$f.eid
colnames(mdd.phenotypes)[grep('f.eid',colnames(mdd.phenotypes))]='FID'
mdd.phenotypes=mdd.phenotypes[,c('FID','IID','mdd_nerves','mdd_smith','mdd_icd','Depressed.Ever')]
colnames(mdd.phenotypes)[ncol(mdd.phenotypes)]='mdd_cidi'

ID.tokeep = read.table('/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/Share/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id',header=F,stringsAsFactors=F)  # unrelated Caucasian sample

mdd.phenotypes.plink=mdd.phenotypes[mdd.phenotypes$IID %in% ID.tokeep$V1,]

write.table(mdd.phenotypes.plink,file='data/ukb_mdd_phenotypes.txt',col.names=T,row.names=F,quote=F,sep=' ')
write.table(mdd.phenotypes.plink[,1:2],file='data/ukb_unrelated_white_ID.txt',col.names=T,row.names=F,quote=F,sep=' ')

# plink-format data for unrelated caucasian sample ------------------------------------------

ukb.genotype='/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_genotype/autosome.qc.maf01.hwe5e-6.geno02.mind02.snps'
target.genotype='/exports/eddie/scratch/xshen33/autosome.qc.maf01.hwe5e-6.geno02.mind02.snps.unrelatedcaucasian'
ID.list='/gpfs/igmmfs01/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/data/ukb_unrelated_white_ID.txt'
system(paste0('plink --bfile ',ukb.genotype,' --keep ',ID.list,' --make-bed --out ',target.genotype))

# bgen-format data matched phenotype files ------------------------------------------
ID.bgen=read.table('/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3/ukb4844_imp_chr22_v3_s487395.sample',header=F,skip=2)
mdd.phenotypes.bgen=merge(ID.bgen,mdd.phenotypes,by.x='V1',by.y='IID',all.x=T)
mdd.phenotypes.bgen=mdd.phenotypes.bgen[,c(1,2,6:9)]
colnames(mdd.phenotypes.bgen)[1:2]=c('FID','IID')
rownames(mdd.phenotypes.bgen)=mdd.phenotypes.bgen$IID
mdd.phenotypes.bgen=mdd.phenotypes.bgen[as.character(ID.bgen$V1),]

write.table(mdd.phenotypes.bgen,file='data/ukb_mdd_phenotypes_bgen.txt',col.names=T,row.names=F,quote=F,sep=' ')