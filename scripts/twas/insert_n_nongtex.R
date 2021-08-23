#!/usr/bin/Rscript
panels<-c("CMC.BRAIN.RNASEQ","CMC.BRAIN.RNASEQ_SPLICING","NTR.BLOOD.RNAARR","YFS.BLOOD.RNAARR")

panels<-data.frame(panel=panels,
                   N=c(452,452,1247,1264))

library(data.table)

for(i in panels$panel){
  pos<-fread(paste0('resources/twas/fusion_data/',i,'/',i,'.pos'))
  pos$PANEL<-i
  pos$N<-panels$N[panels$panel == i]
  pos<-pos[,c('PANEL','WGT','ID','CHR','P0','P1','N'), with=F]
  write.table(pos, paste0('resources/twas/fusion_data/',i,'/',i,'.pos'), quote=F, col.names=T, row.names=F)
}
