#!/usr/bin/Rscript

# Change so the WGT column is used instead of the ID column
# This makes the FOCUS results more distinguishable before panels and splice variants

library(data.table)

# Start with the FUSION panels
panels<-c('Adrenal_Gland','Brain_Amygdala','Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus','Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia','Brain_Substantia_nigra','CMC.BRAIN.RNASEQ','CMC.BRAIN.RNASEQ_SPLICING','NTR.BLOOD.RNAARR','Pituitary','Thyroid','Whole_Blood','YFS.BLOOD.RNAARR')

for(i in panels){
  pos_i<-fread(paste0('resources/twas/fusion_data/',i,'/',i,'.pos'))
  pos_i$ID<-pos_i$WGT
  write.table(pos_i, paste0('resources/twas/fusion_data/',i,'/',i,'.WGT_ID.pos'), col.names=T, row.names=F, quote=F)
}

# Some PsychENCODE SNP-weights don't have any weights. Make a .pos file that does not include these.
# Also create a version with WGT as ID
pos<-read.table('resources/twas/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights.pos', header=T, stringsAsFactors=F)

pos_clean<-NULL
for(i in 1:dim(pos)[1]){
  load(paste0('resources/twas/psychencode_data/PEC_TWAS_weights/',pos$WGT[i]))
  
  print(i)
  if(dim(snps)[1] > 0){
    pos_clean<-rbind(pos_clean,pos[i,])
  }
}

dim(pos)
dim(pos_clean)
# They seem to have resolved this issue now, but worth checking.

pos_clean_WGT<-pos_clean
pos_clean_WGT$ID<-pos_clean_WGT$WGT

write.table(pos_clean, 'resources/twas/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights_noMissing.pos', col.names=T, row.names=F, quote=F)

write.table(pos_clean_WGT, 'resources/twas/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights_noMissing.WGT_ID.pos', col.names=T, row.names=F, quote=F)

# Create FOCUS database
for(panel in panels){
  system(paste0('focus import resources/twas/fusion_data/',panel,'/',panel,'.WGT_ID.pos fusion --tissue ',panel,' --output resources/twas/focus_db/MDD_TWAS'))
}

system('focus import resources/twas/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights_noMissing.WGT_ID.pos fusion --tissue PEC_TWAS_weights --output resources/twas/focus_db/MDD_TWAS')
