#!/usr/bin/Rscript

library(data.table)
twas <- fread("results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt")
twas_sign <- fread("results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_TWSig.txt")

cat(dim(twas)[1], 'features were present in SNP-weight sets\n')
cat(sum(!is.na(twas$TWAS.P)), 'features were imputed in TWAS\n')
cat(length(unique(twas$ID)), 'unique features were present in SNP-weight sets\n')
cat(length(unique(twas$ID[!is.na(twas$TWAS.P)])), 'unique features were imputed in TWAS\n')

twas_sign$BEST.GWAS.P<-2*pnorm(-abs(twas_sign$BEST.GWAS.Z))
cat(sum(twas_sign$BEST.GWAS.P > 5e-8), 'TWSig features are not within GWSig SNPs\n')
cat(sum(twas_sign$BEST.GWAS.P < 5e-8), 'TWSig features are within GWSig SNPs\n')

# Clean the PANEL names of the output df containing results on all tested features
twas$PANEL_clean<-gsub('_',' ',twas$PANEL)
twas$PANEL_clean<-gsub('CMC.BRAIN.RNASEQ','CMC DLPFC',twas$PANEL_clean)
twas$PANEL_clean<-gsub('SPLICING','Splicing',twas$PANEL_clean)
twas$PANEL_clean<-gsub('NTR.BLOOD.RNAARR','NTR Blood',twas$PANEL_clean)
twas$PANEL_clean<-gsub('YFS.BLOOD.RNAARR','YFS Blood',twas$PANEL_clean)
twas$PANEL_clean[!grepl('CMC|NTR|YFS|PsychENCODE', twas$PANEL)]<-paste0('GTEx ',twas$PANEL_clean[!grepl('CMC|NTR|YFS|PsychENCODE', twas$PANEL)])
#to add gtex to each of the snp weights which don't have CMC NTR or YFS in front
twas$PANEL_clean<-gsub('Brain', '', twas$PANEL_clean)
twas$PANEL_clean <- gsub('Anterior cingulate cortex', 'ACC', twas$PANEL_clean)
twas$PANEL_clean <- gsub('basal ganglia', '', twas$PANEL_clean)
twas$PANEL_clean <- gsub('BA9', '', twas$PANEL_clean)
twas$PANEL_clean <- gsub('BA24', '', twas$PANEL_clean)
twas$PANEL_clean <- gsub('  ', ' ', twas$PANEL_clean)

# Create a table showing the number of features tested for each SNP-weight set
tab_ob<-table(twas$PANEL_clean)
panel_tab<-data.frame(tab_ob)
names(panel_tab)<-c('PANEL','N_feat')

tab_imp_ob<-table(twas[!is.na(twas$TWAS.P),]$PANEL_clean)
panel_imp_tab<-data.frame(tab_imp_ob)
names(panel_imp_tab)<-c('PANEL','N_feat')

panel_tab_all<-merge(panel_tab, panel_imp_tab, by='PANEL')
names(panel_tab_all)<-c('PANEL','N_feat','N_feat_imp')

write.csv(panel_tab_all, file = "results/twas/twas_results/PGC_MDD3_twas_panel_N.csv", row.names = F)

# Shorten panel name to plot easily
twas$PANEL_clean_short<-substr(twas$PANEL_clean, start = 1, stop = 12)  #start the name at the first character and stop at the 25th
twas$PANEL_clean_short[nchar(twas$PANEL_clean) > 12]<-paste0(twas$PANEL_clean_short[nchar(twas$PANEL_clean) > 12], "...")

###
# Deal with missingness and subset for the relevant cols only
###

##TWAS df
#exclude missings
twas<-twas[!is.na(twas$TWAS.Z),]
twas<-twas[!is.na(twas$TWAS.P),]

#subset columns needed 
twas_sub <- twas[,c('FILE', 'ID','PANEL', 'PANEL_clean_short','PANEL_clean','CHR','P0', 'P1', 'TWAS.Z', 'TWAS.P', 'COLOC.PP0', 'COLOC.PP1', 'COLOC.PP2', 'COLOC.PP3', 'COLOC.PP4')]

###
# Update positions
###

# Rationale: the positions in the output files created by FUSION are rounded, thus not completely accurate. 
# Therefore, we need to update the positions (P0 and P1) based on the pos files in Rosalind. 

weights<-c('Adrenal_Gland', 'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Brain_Substantia_nigra', 'CMC.BRAIN.RNASEQ', 'CMC.BRAIN.RNASEQ_SPLICING', 'NTR.BLOOD.RNAARR', 'Pituitary', 'Thyroid', 'Whole_Blood', 'YFS.BLOOD.RNAARR')

#Get all pos files within the SNP-weight sets and bind them 
FUSION_pos<-NULL
for(i in weights){
  FUSION_pos_temp<-read.table(paste('resources/twas/fusion_data/',i, '/',i, '.pos',sep=''), header=T, stringsAsFactors=F)   #repeating i twice with / in the middle is to get one folder further
  FUSION_pos<-rbind(FUSION_pos, FUSION_pos_temp)
}

PsychENCODE_pos<-read.table('resources/twas/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights.pos', header=T, stringsAsFactors=F)

# Combine pos files
FUSION_pos<-rbind(FUSION_pos, PsychENCODE_pos)

write.table(FUSION_pos,'results/twas/twas_results/PGC_MDD3_twas.pos', col.names=T, row.names=F, quote=F)

###
# Merge the pos file with the twas_sub
###

#the pos file and the output file do not have the same columns with the same information. We therefore need to slightly modify the TWAS columns 
twas_sub$tmp<-gsub('.*resources/twas/fusion_data/','',twas_sub$FILE)
twas_sub$tmp<-gsub('.*resources/twas/psychencode_data/','',twas_sub$tmp)
#to delete the full pathway of the file and just keep the important information 
twas_sub$PANEL<-sub('/.*','', twas_sub$tmp)
twas_sub$Feature<-gsub('.*/','',twas_sub$tmp)
twas_sub$WGT<-paste0(twas_sub$PANEL, '/', twas_sub$Feature)
twas_sub$PANEL<-NULL
twas_sub$tmp<-NULL
twas_sub$Feature<-NULL

#merge
twas_sub_correct <- merge(twas_sub, FUSION_pos, by="WGT")

#clean
twas_sub_correct$ID.y<-NULL
names(twas_sub_correct)[3]<-'ID'   #to change the name  of IDx to ID

###
# Clean output files for future scripts
###

names(twas_sub_correct)[names(twas_sub_correct) == 'CHR.x']<- "CHR"
names(twas_sub_correct)[names(twas_sub_correct) == 'P0.y']<- "P0"
names(twas_sub_correct)[names(twas_sub_correct) == 'P1.y']<- "P1"
twas_sub_correct$CHR.y <- NULL
twas_sub_correct$P0.x <- NULL
twas_sub_correct$P1.x <- NULL

twas_sub_correct$N <- NULL

#save 
write.table(twas_sub_correct, file = "results/twas/twas_results/PGC_MDD3_twas_AllTissues_CLEAN.txt", sep = " ", col.names = T, row.names = F)
write.table(twas_sub_correct[twas_sub_correct$TWAS.P < 1.368572e-06,], file = "results/twas/twas_results/PGC_MDD3_twas_AllTissues_TWSig_CLEAN.txt", sep = " ", col.names = T, row.names=F)
