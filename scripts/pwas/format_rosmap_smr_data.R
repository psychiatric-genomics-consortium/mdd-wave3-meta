#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--rosmap", action="store", default=NA, type='character',
              help="Path to ROSMAP data [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

pQTL<-fread(opt$rosmap)

# Convert to matrix QTL format
pQTL$T_stat<-pQTL$BETA/pQTL$SE
pQTL$FDR<-p.adjust(pQTL$P, method='fdr')

pQTL_mat<-pQTL[,c('SNP','Protein_UniProt','Beta','T_stat','P','FDR')]
names(pQTL_mat)<-c('SNP','gene','beta','t-stat','p-value','FDR')

dir.create('resources/data/rosmap_smr/')

fwrite(pQTL_mat, 'resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt', sep=' ', quote=F, na='NA')

# Convert to BESD format for SMR
system(paste0('resources/twas/smr/smr_Linux --eqtl-summary resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt --matrix-eqtl-format --make-besd --out resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd'))

# Update the .esi files
pQTL<-fread(opt$rosmap)
esi<-fread('resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd.esi')

pQTL_unique_snp<-pQTL[!duplicated(pQTL$SNP),]
esi_new<-merge(esi, pQTL_unique_snp, by.x='V2', by.y='SNP')
esi_new<-esi_new[,c('Chr','V2','V3','BP','A1','A2'),with=F]
names(esi_new)<-c('CHR','SNP','POS','BP','A1','A2')

frq<-NULL
for(i in 1:22){
  tmp<-fread(paste0('/users/k1806347/brc_scratch/Software/mdd-meta/resources/twas/1kg/EUR/EUR_phase3.MAF_001.chr',i,'.frq'))
  tmp<-tmp[tmp$SNP %in% esi_new$SNP[esi_new$CHR == i],]
  frq<-rbind(frq, tmp)
}

esi_new_match<-merge(esi_new, frq[,c('SNP','A1','A2','MAF'),with=F], by=c('SNP','A1','A2'))
esi_new_swap<-merge(esi_new, frq[,c('SNP','A1','A2','MAF'),with=F], by.x=c('SNP','A1','A2'), by.y=c('SNP','A2','A1'))
esi_new_swap$MAF<-1-esi_new_swap$MAF

esi_new_unmatch<-merge(esi_new, frq[,c('SNP','MAF'),with=F], by='SNP',all.x=T)
esi_new_unmatch<-esi_new_unmatch[!(esi_new_unmatch$SNP %in% esi_new_match$SNP) & !(esi_new_unmatch$SNP %in% esi_new_swap$SNP)]
# Insert rs12186596 MAF as 0.3557 (TWINSUK estimate)
esi_new_unmatch$MAF<-0.3557

# Only one SNP isn't in the reference
esi_new_freq<-do.call(rbind, list(esi_new_match,esi_new_swap,esi_new_unmatch))
esi_new_freq<-esi_new_freq[,c('CHR','SNP','POS','BP','A1','A2','MAF'),with=F]

fwrite(esi_new_freq, 'resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd.esi_update', sep=' ', quote=F, na='NA')

###
# Update the .epi file
###

pQTL_unique_gene<-pQTL[!duplicated(pQTL$Protein_UniProt),]
epi<-fread('resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd.epi')

epi_new<-merge(epi, pQTL_unique_gene[,c('Protein_Chr','Protein_UniProt','Protein_BP_Start','Protein_GeneSymbol'), with=F], by.x='V2', by.y='Protein_UniProt')

epi_new$Dir<-'+'
epi_new$POS<-0

epi_new<-epi_new[,c('Protein_Chr','V2','POS','Protein_BP_Start','Protein_GeneSymbol','Dir'),with=F]

fwrite(epi_new, 'resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd.epi_update', sep=' ', quote=F, na='NA')

# Update the esi
system(paste0('resources/twas/smr/smr_Linux --beqtl-summary resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd --update-esi resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd.esi_update'))

# Update the epi
system(paste0('resources/twas/smr/smr_Linux --beqtl-summary resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd --update-epi resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd.epi_update'))


