#!/usr/bin/Rscript

library(data.table)
res<-fread('results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt')
Brain_res<-res[grepl('Brain|BRAIN|PsychENCODE', res$PANEL),]
HPA_res<-res[grepl('Adrenal|Pituitary|Hypothalamus', res$PANEL),]
HPT_res<-res[grepl('Thyroid|Pituitary|Hypothalamus', res$PANEL),]
BLOOD_res<-res[grepl('BLOOD', res$PANEL),]

write.table(Brain_res, 'results/twas/twas_results/PGC_MDD3_twas_BRAIN_GW.txt', row.names=F, col.names=T, quote=F)
write.table(HPA_res, 'results/twas/twas_results/PGC_MDD3_twas_HPA_GW.txt', row.names=F, col.names=T, quote=F)
write.table(HPT_res, 'results/twas/twas_results/PGC_MDD3_twas_HPT_GW.txt', row.names=F, col.names=T, quote=F)
write.table(BLOOD_res, 'results/twas/twas_results/PGC_MDD3_twas_BLOOD_GW.txt', row.names=F, col.names=T, quote=F)