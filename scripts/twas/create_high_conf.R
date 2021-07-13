#!/usr/bin/Rscript

library(data.table)

focus<-fread('results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.results.csv')
fusion<-fread('results/twas/conditional/PGC_MDD3_TWAS_Conditional_table_novelty.csv')

fusion_focus<-merge(fusion, focus[,c('mol_name','tissue','twas_z','pip','region'),with=F], by.x=c('WGT'), by.y=c('mol_name'), all.x=T)
fusion_focus<-fusion_focus[,c('WGT','CHR','P0','P1','PANEL_clean','ID','TWAS.Z','TWAS.P','Novel','COLOC.PP4','Colocalised','pip','region'),with=F]
names(fusion_focus)<-c('WGT','CHR','P0','P1','SNP-weight Set','ID','TWAS.Z','TWAS.P','Novel','COLOC.PP4','Colocalised','FOCUS_pip','FOCUS_region')
fusion_focus<-fusion_focus[order(fusion_focus$CHR, fusion_focus$P0),]
fusion_focus$Location<-paste0('chr',fusion_focus$CHR,':',fusion_focus$P0,'-',fusion_focus$P1)   

# Remove the MHC region
fusion_focus_noMHC<-fusion_focus[!(fusion_focus$CHR == 6 & fusion_focus$P1 > 26e6 & fusion_focus$P0 < 34e6),]

# Subset those which are high confidence
fusion_focus_highConf<-fusion_focus_noMHC[fusion_focus_noMHC$Colocalised == T & fusion_focus_noMHC$FOCUS_pip > 0.5 & fusion_focus_noMHC$TWAS.P < 3.685926e-08,]

write.csv(fusion_focus_highConf,'results/twas/PGC3_MDD_TWAS_HighConf_results.csv', row.names=F, quote=F)
