#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--out", action="store", default=NA, type='character',
              help="Name of output file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(biomaRt)

weights<-c("Adrenal_Gland","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Substantia_nigra","CMC.BRAIN.RNASEQ","CMC.BRAIN.RNASEQ_SPLICING","NTR.BLOOD.RNAARR","Pituitary","Thyroid","Whole_Blood","YFS.BLOOD.RNAARR")

# Write out this list of SNP-weights as this might be useful elsewhere
write.table(c(weights,'PsychENCODE'), '/users/k1806347/brc_scratch/Software/mdd-meta/results/twas/list_of_weights.txt', col.names=F, row.names=F, quote=F) 

all<-NULL

###
# Read in and format PsychENCODE twas results
###

psychENCODE<-NULL
for(i in 1:22){
  if(i == 6){
    tmp1<-fread(paste0('results/twas/psychencode/PGC_MDD3_twas_psychencode_chr',i))
    tmp2<-fread(paste0('results/twas/psychencode/PGC_MDD3_twas_psychencode_chr',i,'.MHC'))
    tmp<-rbind(tmp1,tmp2)
  } else {
    tmp<-fread(paste0('results/twas/psychencode/PGC_MDD3_twas_psychencode_chr',i))
  }
  
  psychENCODE<-rbind(psychENCODE, tmp)
}

col_order<-names(psychENCODE)

psychENCODE$PANEL<-as.character(psychENCODE$PANEL)
psychENCODE$PANEL<-'PsychENCODE'

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = ensembl)

psychENCODE<-merge(psychENCODE, Genes, by.x='ID', by.y='ensembl_gene_id')
psychENCODE$ID<-psychENCODE$external_gene_name
psychENCODE$external_gene_name<-NULL

psychENCODE<-psychENCODE[,col_order, with=F]

write.table(psychENCODE, paste0(opt$out,'PGC_MDD3_twas_PsychENCODE_GW.txt'), col.names=T, row.names=F, quote=F)

all<-rbind(all, psychENCODE)

###
# Read in and format FUSION twas results
###

for(weight in weights){
  FUSION<-NULL
  for(i in 1:22){
    if(i == 6){
      tmp1<-fread(paste0('results/twas/PGC_MDD3_twas_',weight,'_chr',i))
      tmp2<-fread(paste0('results/twas/PGC_MDD3_twas_',weight,'_chr',i,'.MHC'))
      tmp<-rbind(tmp1,tmp2)
    } else {
      tmp<-fread(paste0('results/twas/PGC_MDD3_twas_',weight,'_chr',i))
    }
    
    FUSION<-rbind(FUSION, tmp)
  }
  write.table(FUSION, paste0(opt$out,'PGC_MDD3_twas_',weight,'_GW.txt'), col.names=T, row.names=F, quote=F)
  all<-rbind(all, FUSION)
}

###
# Write out combined results
###

write.table(all, paste0(opt$out,'PGC_MDD3_twas_AllTissues_GW.txt'), row.names=F, col.names=T, quote=F)

# Write out transcriptome-wide significant results
write.table(all[which(all$TWAS.P < 1.368572e-06),], paste0(opt$out,'PGC_MDD3_twas_AllTissues_GW_TWSig.txt'), row.names=F, col.names=T, quote=F)

cat('N TWAS-sig features:',dim(all[which(all$TWAS.P < 1.368572e-06),]),'\n',sep='')
cat('N TWAS-sig unique genes:',length(unique(all[which(all$TWAS.P < 1.368572e-06),]$ID)),'\n',sep='')
