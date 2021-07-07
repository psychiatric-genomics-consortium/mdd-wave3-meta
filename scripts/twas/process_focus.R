#!/usr/bin/Rscript

library(data.table)

fusion <- fread("results/twas/twas_results/PGC_MDD3_twas_AllTissues_TWSig_CLEAN.txt")

focus.files<-list.files(path='results/twas/focus/p_1e-4/', pattern=glob2rx("PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr*.focus.tsv"))
length(focus.files)
focus<-NULL
for(i in focus.files){
focus<-rbind(focus,fread(paste0('results/twas/focus/p_1e-4/',i)))
}

# Update the feature IDs with gene names
pos<-fread('results/twas/twas_results/PGC_MDD3_twas.pos')
focus<-merge(focus, pos[,c('WGT','ID')], by.x='mol_name', by.y='WGT',all.x=T)
focus<-focus[order(focus$chrom, focus$region, 1-focus$pip),]

# I noticed a bug in the output where features that should be in the 90% credible set are not
focus_bug<-NULL
for(i in unique(focus$region)){
    focus_temp<-focus[focus$region == i,]
    if(sum(focus_temp$in_cred_set) == 0 & max(focus_temp$pip) != focus_temp$pip[focus_temp$ens_gene_id == 'NULL.MODEL']){
        focus_bug<-rbind(focus_bug, focus_temp)
    }
	focus$in_cred_set[focus$ID == focus_temp$ID[1] & focus$tissue == focus_temp$tissue[1] & focus$region == focus_temp$region[1]] <- 1
}

# Update tissue for psychencode features
focus_psychencode<-focus[focus$tissue == 'pec_twas_weights',]
focus_fusion<-focus[focus$tissue != 'pec_twas_weights',]
focus_psychencode$tissue<-'psychencode'

# Update PsychENCODE gene IDs from ensembl to gene names
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = ensembl)

focus_psychencode<-merge(focus_psychencode, Genes, by.x='ID', by.y='ensembl_gene_id')
focus_psychencode$ID<-focus_psychencode$external_gene_name
focus_psychencode$external_gene_name<-NULL
focus_psychencode<-focus_psychencode[,names(focus),with=F]
focus<-rbind(focus_fusion,focus_psychencode)

fusion_focus<-merge(fusion, focus[,c('mol_name','tissue','twas_z','pip','in_cred_set','region'),with=F], by.x=c('WGT'), by.y=c('mol_name'), all.x=T)
fusion_focus<-fusion_focus[,c('WGT','CHR','P0','P1','PANEL_clean_short','ID','TWAS.Z','TWAS.P','twas_z','in_cred_set','pip','region'),with=F]
names(fusion_focus)<-c('WGT','CHR','P0','P1','SNP-weight Set','ID','TWAS.Z','TWAS.P','FOCUS_twas_z','FOCUS_in_cred_set','FOCUS_pip','FOCUS_region')
fusion_focus<-fusion_focus[order(fusion_focus$CHR, fusion_focus$P0),]
fusion_focus$Location<-paste0('chr',fusion_focus$CHR,':',fusion_focus$P0,'-',fusion_focus$P1)   
fusion_focus<-fusion_focus[,c('Location','SNP-weight Set','ID','TWAS.Z','TWAS.P','FOCUS_twas_z','FOCUS_in_cred_set','FOCUS_pip','FOCUS_region'),with=F]

write.csv(fusion_focus,'results/twas/focus/p_1e-4/PGC_MDD3_TWAS.TWSig.FOCUS.results.csv', row.names=F, quote=F)
write.csv(focus,'results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.results.csv', row.names=F, quote=F)