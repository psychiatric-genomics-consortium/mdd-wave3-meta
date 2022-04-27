#!/usr/bin/Rscript

library(data.table)
library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = ensembl)

Genes<-Genes[!duplicated(Genes$external_gene_name),]

# Read in the rosmap smr results
smr_rosmap_files<-list.files(path='results/pwas/rosmap_smr/', pattern='rosmap_smr_res_chr')
smr_rosmap_files<-smr_rosmap_files[grepl('.smr$', smr_rosmap_files)]

smr_rosmap<-NULL
for(i in smr_rosmap_files){
  smr_rosmap<-rbind(smr_rosmap, fread(paste0('results/pwas/rosmap_smr/',i)))
}

# Split rows containing a string of gene names into seperate rows
smr_rosmap_one_gene<-smr_rosmap[!grepl(';', smr_rosmap$Gene),]
smr_rosmap_mult_gene<-smr_rosmap[grepl(';', smr_rosmap$Gene),]

smr_rosmap_mult_gene_split<-NULL
for(i in 1:nrow(smr_rosmap_mult_gene)){
  ids<-unlist(strsplit(smr_rosmap_mult_gene$Gene[i], ';'))
  for(j in ids){
    tmp<-smr_rosmap_mult_gene[i,]
    tmp$Gene<-j
    smr_rosmap_mult_gene_split<-rbind(smr_rosmap_mult_gene_split, tmp)
  }
}

smr_rosmap<-rbind(smr_rosmap_one_gene, smr_rosmap_mult_gene_split)

smr_rosmap<-smr_rosmap[!is.na(smr_rosmap$Gene),]
smr_rosmap<-smr_rosmap[!duplicated(smr_rosmap$Gene),]
smr_rosmap$external_gene_name<-smr_rosmap$Gene

smr_rosmap<-merge(smr_rosmap,Genes, by='external_gene_name', all.x=T)
smr_rosmap<-smr_rosmap[!duplicated(smr_rosmap$ensembl_gene_id),]

fwrite(smr_rosmap, 'results/pwas/rosmap_smr/rosmap_smr_res_GW.txt.gz', quote=F, sep=' ', na='NA')

