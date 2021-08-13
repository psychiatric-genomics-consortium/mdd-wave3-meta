#!/usr/bin/Rscript

# Create a list of ensemble IDs
IDs<-list.files('resources/twas/psychencode_data/SNP-weights/PEC_TWAS_weights')
IDs<-IDs[grepl('.wgt.RDat', IDs)]
IDs<-gsub('.wgt.RDat','',IDs)

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), mart = ensembl)
Genes<-Genes[(Genes$ensembl_gene_id %in% IDs),]
Genes$chromosome_name<-as.numeric(Genes$chromosome_name)
Genes<-Genes[,c("chromosome_name","start_position","end_position","ensembl_gene_id")]
Genes<-Genes[order(Genes$chromosome_name, Genes$start_position),]

write.table(Genes, 'resources/twas/psychencode_data/PEC_TWAS_weights.coord', col.names=T, row.names=F, quote=F)

system('Rscript resources/twas/Calculating-FUSION-TWAS-weights-pipeline/OP_packaging_fusion_weights.R --RDat_dir resources/twas/psychencode_data/SNP-weights/PEC_TWAS_weights --coordinate_file resources/twas/psychencode_data/PEC_TWAS_weights.coord --output_name PEC_TWAS_weights --output_dir resources/twas/psychencode_data/PEC_TWAS_weights')

pos<-read.table('resources/twas/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights.pos', header=T)
pos$N<-1321
write.table(pos, 'resources/twas/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights.pos', col.names=T, row.names=F, quote=F)
