#!/usr/bin/Rscript
library(data.table)

tissue<-c('Basalganglia','Cerebellum','Cortex','Hippocampus','Spinalcord')
chr<-1:22

for(tissue_i in tissue){
  for(chr_i in chr){
    esi<-fread(paste0('resources/twas/MetaBrain/',tissue_i,'/2020-05-26-',tissue_i,'-EUR-',chr_i,'-SMR-besd.esi'))
    esi$V2_new<-esi$V2
    esi$V2_new<-gsub(':.*','',gsub('.*rs','rs',esi$V2_new))
    esi$V2_new[esi$V2_new == 'rs']<-'nors'
    esi$V2<-esi$V2_new
    esi$V2_new<-NULL
    fwrite(esi, paste0('resources/twas/MetaBrain/',tissue_i,'/2020-05-26-',tissue_i,'-EUR-',chr_i,'-SMR-besd.esi'), col.names=F, row.names=F, quote=F, sep='\t', na='NA')
  }
}
