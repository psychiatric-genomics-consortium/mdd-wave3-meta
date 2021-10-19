#!/usr/bin/Rscript
library(data.table)

#####
# eQTLGen
#####
smr_res<-NULL
for(i in 1:22){
  smr_res<-rbind(smr_res,fread(paste0('results/twas/eqtlgen_smr/eqtlgen_smr_res_chr',i,'.smr')))
}

# Insert gene names
# The Ensembl version 71 (used by eQTLGen) is not available on biomaRt so 89 genes do not have IDs
# To get a perfect match, download the significant eQTL table to link Ensembl IDs to gene names
setwd('resources/twas')
system(paste0('wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz'))
eqtlres<-fread('2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')
eqtlres<-eqtlres[!duplicated(eqtlres$Gene),]
eqtlres<-eqtlres[,c('Gene','GeneSymbol'), with=F]
system(paste0('rm 2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz'))
setwd('../..')

smr_res<-merge(smr_res, eqtlres, by.x='probeID', by.y='Gene', all.x=T)

# Calculate Bonferroni corrected p-values
smr_res$p_SMR_bonf<-p.adjust(smr_res$p_SMR, method="bonferroni")

# Save the full SMR results
write.csv(smr_res, 'results/twas/eqtlgen_smr/eqtlgen_smr_res_GW_withIDs.csv', row.names=F, quote=F)

######
# MetaBrain
######

tissue<-c("Basalganglia","Cerebellum","Cortex","Hippocampus","Spinalcord")

for(tissue_i in tissue){
  smr_res<-NULL
  for(chr_i in 1:22){
    smr_res<-rbind(smr_res,fread(paste0('results/twas/metabrain_smr/metabrain_',tissue_i,'_smr_res_chr',chr_i,'.smr')))
  }
  
  # Again, some Ensembl IDs aren't coming up in BioMart, so lets download the eQTL data on the MetaBrain website which contain gene symbols.
  meta_brain_tmp<-NULL
  for(chr_i in 1:22){
    setwd('resources/twas')
    system(paste0('wget https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-',tissue_i,'-EUR/2020-05-26-',tissue_i,'-EUR-',chr_i,'-biogenformat.txt.gz'))
    setwd('../..')
    tmp<-fread(paste0('resources/twas/2020-05-26-',tissue_i,'-EUR-',chr_i,'-biogenformat.txt.gz'), select=c('ProbeName','HGNCName'))
    tmp<-tmp[!duplicated(tmp),]
    meta_brain_tmp<-rbind(meta_brain_tmp, tmp)
    system(paste0('rm resources/twas/2020-05-26-',tissue_i,'-EUR-',chr_i,'-biogenformat.txt.gz'))
  }
  
  meta_brain_tmp<-meta_brain_tmp[!duplicated(meta_brain_tmp),]
  
  smr_res<-merge(smr_res, meta_brain_tmp, by.x='Gene', by.y='ProbeName', all.x=T)
  
  # Calculate Bonferroni corrected p-values
  smr_res$p_SMR_bonf<-p.adjust(smr_res$p_SMR, method="bonferroni")
  
  # Save the full SMR results
  write.csv(smr_res, paste0('results/twas/metabrain_smr/metabrain_',tissue_i,'_smr_res_GW_withIDs.csv'), row.names=F, quote=F)
  
}
