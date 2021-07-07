#!/usr/bin/Rscript

# Read in the report files
library(data.table)

# Read in the clean TWAS results
twas_sign <- fread("results/twas/twas_results/PGC_MDD3_twas_AllTissues_TWSig_CLEAN.txt")
twas_sign$PANEL_clean<-gsub(' $','',twas_sign$PANEL_clean)

# Read in all jointly significant associations
temp = list.files(path='results/twas/conditional/',pattern=glob2rx("*chr*.report"))
report<-do.call(rbind, lapply(temp, function(x) read.table(paste0('results/twas/conditional/',x), header=T,stringsAsFactors=F)))
report$JOINT.ID<-NA
report$MARGIN.ID<-NA
report$JOINT.N<-NA
report$MARGIN.N<-NA
report$loc<-gsub('.*loc_','',report$FILE)
joint_res<-NULL
margin_res<-NULL

# Insert names of jointly significant genes
for(i in unique(report$CHR)){
  joint_i<-read.table(paste0('results/twas/conditional/PGC_MDD3_TWAS_conditional_chr',i,'.joint_included.dat'), header=T,stringsAsFactors=F)
  margin_i<-read.table(paste0('results/twas/conditional/PGC_MDD3_TWAS_conditional_chr',i,'.joint_dropped.dat'), header=T,stringsAsFactors=F)
  
  joint_i$path<-gsub('/[^/]+$','',joint_i$FILE)
    joint_i$path<-gsub('/[^/]+$','',joint_i$path)
    joint_i$WGT<-NA
    for(j in 1:dim(joint_i)[1]){
      joint_i$WGT[j]<-gsub(paste0(joint_i$path[j],'/'),'',joint_i$FILE[j])
  }

  if(dim(margin_i)[1] > 0){
    margin_i$path<-gsub('/[^/]+$','',margin_i$FILE)
    margin_i$path<-gsub('/[^/]+$','',margin_i$path)
    margin_i$WGT<-NA
    for(j in 1:dim(margin_i)[1]){
      margin_i$WGT[j]<-gsub(paste0(margin_i$path[j],'/'),'',margin_i$FILE[j])
    }
  }

  temp = list.files(path="results/twas/conditional/", pattern=glob2rx(paste0("*chr",i,".loc*.genes")))

  for(k in 1:length(temp)){
    loc_k<-read.table(paste0('results/twas/conditional/PGC_MDD3_TWAS_conditional_chr',i,'.loc_',k,'.genes'), header=T, stringsAsFactors=F)
    
    loc_k$path<-gsub('/[^/]+$','',loc_k$FILE)
    loc_k$path<-gsub('/[^/]+$','',loc_k$path)
    loc_k$WGT<-NA
    for(j in 1:dim(loc_k)[1]){
      loc_k$WGT[j]<-gsub(paste0(loc_k$path[j],'/'),'',loc_k$FILE[j])
    }

    loc_k$P0<-NULL
    loc_k$P1<-NULL
    
    loc_k<-merge(loc_k, twas_sign[,c('P0','P1','WGT','PANEL_clean')], by='WGT')
    
    loc_k_joint<-loc_k[(loc_k$WGT %in% joint_i$WGT),]
    joint_res<-rbind(joint_res,loc_k_joint)
        
    if(dim(margin_i)[1] > 0){
      loc_k_margin<-loc_k[(loc_k$WGT %in% margin_i$WGT),]
      margin_res<-rbind(margin_res,loc_k_margin)
    } else {
      loc_k_margin<-data.frame(ID=NULL)
    }
    
    g_list<-NULL
    for(g in unique(loc_k_joint$ID)){
      g_list<-c(g_list,paste0(g, " (",paste(loc_k_joint$PANEL_clean[loc_k_joint$ID == g], collapse=', '),")"))
    }
    report[report$CHR == i & report$loc == k,]$JOINT.ID<-paste(g_list,collapse=', ')

    if(dim(loc_k_margin)[1] > 0){
      g_list<-NULL
      for(g in unique(loc_k_margin$ID)){
        g_list<-c(g_list,paste0(g, " (",paste(unique(loc_k_margin$PANEL_clean[loc_k_margin$ID == g]), collapse=', '),")"))
      }
      report[report$CHR == i & report$loc == k,]$MARGIN.ID<-paste(g_list,collapse=', ')
    } else {
      report[report$CHR == i & report$loc == k,]$MARGIN.ID<-'-'
    }
    
    report[report$CHR == i & report$loc == k,]$JOINT.N<-dim(loc_k_joint)[1]
    report[report$CHR == i & report$loc == k,]$MARGIN.N<-dim(loc_k_margin)[1]
  }
}

report$LOCUS<-paste0(report$CHR,':',report$P0,':',report$P1)
report$BP<-paste0(report$P0,'-',report$P1)
report$VAR.EXP<-paste0(report$VAR.EXP*100,'%')

report<-report[,c('CHR','P0','P1','BP','LOCUS',"JOINT.N",'MARGIN.N','BEST.TWAS.P','BEST.SNP.P','VAR.EXP','JOINT.ID','MARGIN.ID')]

report<-report[order(report$CHR, report$P0),]

# Save full conditional results table
write.csv(report[,c("CHR","BP","JOINT.ID","MARGIN.ID","BEST.TWAS.P","BEST.SNP.P","VAR.EXP")],'results/twas/conditional/PGC_MDD3_TWAS_Conditional_table_full.csv', row.names=F, quote=T)

# Save brief conditional results table
write.csv(report[,c('CHR','BP','JOINT.ID','MARGIN.N','BEST.TWAS.P','BEST.SNP.P','VAR.EXP')],'results/twas/conditional/PGC_MDD3_TWAS_Conditional_table_brief.csv', row.names=F, quote=T)

# Combine gene results for marginal and joint genes
joint_res$Type<-'Joint'
margin_res$Type<-'Marginal'

gene_res<-rbind(joint_res, margin_res)

# Check number of independent associations
cat(dim(joint_res)[1], 'independent associations\n.')

# Check number of independent associations without genome-wide significant snp
cat(dim(joint_res[2*pnorm(-abs(joint_res$BEST.GWAS.Z)) > 5e-8,])[1], 'independent associations without genome-wide significant snp\n.')

# Check number of independent associations with genome-wide significant snp but an r2 with predicted expression <0.1
cat(dim(joint_res[2*pnorm(-abs(joint_res$BEST.GWAS.Z)) < 5e-8 & joint_res$TOP.SNP.COR^2 < 0.1,])[1], 'independent associations with genome-wide significant snp but an r2 with predicted expression <0.1\n.')

# Check number of independent novel associations
cat(dim(joint_res[(2*pnorm(-abs(joint_res$BEST.GWAS.Z)) < 5e-8 & joint_res$TOP.SNP.COR^2 < 0.1) | 2*pnorm(-abs(joint_res$BEST.GWAS.Z)) > 5e-8,])[1], 'independent novel associations\n.')

# Check number of novel associations
cat(dim(gene_res[(2*pnorm(-abs(gene_res$BEST.GWAS.Z)) < 5e-8 & gene_res$TOP.SNP.COR^2 < 0.1) | 2*pnorm(-abs(gene_res$BEST.GWAS.Z)) > 5e-8,])[1], 'novel associations\n.')

gene_res$Novel<-'No'
gene_res$Novel[(2*pnorm(-abs(gene_res$BEST.GWAS.Z)) < 5e-8 & gene_res$TOP.SNP.COR^2 < 0.1) | 2*pnorm(-abs(gene_res$BEST.GWAS.Z)) > 5e-8]<-'Yes'

gene_res$BP<-paste0(gene_res$P0,'-',gene_res$P1)
gene_res$BEST.GWAS.P<-2*pnorm(-abs(gene_res$BEST.GWAS.Z))

gene_res<-gene_res[order(gene_res$CHR, gene_res$P0),]

gene_res$Colocalised<-F
gene_res$Colocalised[gene_res$COLOC.PP4 >0.8]<-T

# Check number of independent novel associations which colocalise for joint genes
joint_res$Colocalised<-F
joint_res$Colocalised[joint_res$COLOC.PP4 >0.8]<-T

cat(dim(joint_res[(2*pnorm(-abs(joint_res$BEST.GWAS.Z)) < 5e-8 & joint_res$TOP.SNP.COR^2 < 0.1 & joint_res$Colocalised == T) | (2*pnorm(-abs(joint_res$BEST.GWAS.Z)) > 5e-8 & joint_res$Colocalised == T),])[1], 'independent novel associations colocalise for joint genes\n.')

# Check number of novel associations which colocalise
cat(dim(gene_res[(2*pnorm(-abs(gene_res$BEST.GWAS.Z)) < 5e-8 & gene_res$TOP.SNP.COR^2 < 0.1 & gene_res$Colocalised == T) | (2*pnorm(-abs(gene_res$BEST.GWAS.Z)) > 5e-8 & gene_res$Colocalised == T),])[1], 'novel associations colocalise\n.')

gene_res<-gene_res[,c('CHR','BP','P0','P1','ID','PANEL_clean','WGT','TWAS.Z','TWAS.P','BEST.GWAS.P','TOP.SNP.COR','Type','Novel','COLOC.PP3','COLOC.PP4','Colocalised')]

# Save table showing whether gene associations are novel
write.csv(gene_res,'results/twas/conditional/PGC_MDD3_TWAS_Conditional_table_novelty.csv', row.names=F, quote=T)