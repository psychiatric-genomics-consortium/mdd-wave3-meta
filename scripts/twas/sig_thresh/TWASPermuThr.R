#!/usr/bin/Rscript
###############################################################
# Estimate transcriptome-wide significance threshold for PGC3 MDD TWAS
###############################################################

suppressMessages(library("optparse"))

option_list = list(
  make_option("--output", action="store", default=NA, type='character',
              help="output name [required]"),
  make_option("--nperm", action="store", default=100, type='numeric',
              help="number of permutations [required]"),
  make_option("--weights", action="store", default=NA, type='character',
              help="path to file containing a list of SNP-weight sets to be included [required]"),
  make_option("--combinations", action="store", default=NA, type='character',
              help="list of weight group combinations to include [optional]"),
  make_option("--ncore", action="store", default=1, type='character',
              help="number of cores for parallel computing [required]"),
  make_option("--seed", action="store", default=1, type='numeric',
              help="Seed number [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Open library for parallel computing
library(foreach)
library(doMC) 

# Set up parallel environment, specifying the nimber of cores to use
registerDoMC(opt$ncore)

# Load libraries for quickly reading files and regression
library(RcppEigen)
library(data.table)

# Read in list of SNP-weight sets to be included
groups<-fread(opt$weights)
sets<-groups$PANEL

# Read in combinations
if(!is.na(opt$combinations)){
  combinations<-readLines(opt$combinations)
}

# Read in predicted expression levels for all features in the FUSION 1KG reference
if(sum(grepl('PsychENCODE', sets)) == 0){
  GeneX_all<-fread('/mnt/lustre/groups/biomarkers-brc-mh/TWAS_resource/FUSION/Predicted_expression/FUSION_1KG/FUSION_1KG_Expr_AllSets.csv.gz')
} else {
  GeneX_FUSION<-fread('/mnt/lustre/groups/biomarkers-brc-mh/TWAS_resource/FUSION/Predicted_expression/FUSION_1KG/FUSION_1KG_Expr_AllSets.csv.gz')
  GeneX_PsychENCODE<-fread('/scratch/groups/biomarkers-brc-mh/TWAS_resource/PsychEncode/Predicted_expression/FeaturePredictions.csv.gz')
  
  GeneX_all<-merge(GeneX_FUSION, GeneX_PsychENCODE, by=c('FID','IID'))
}

# Read in WGT column from pos files for all SNP-weight sets
pos_all<-NULL
for(i in sets){
  if(i == 'PsychENCODE'){
    pos_temp<-data.frame(WGT=fread('/scratch/groups/biomarkers-brc-mh/TWAS_resource/PsychEncode/PEC_TWAS_weights/PEC_TWAS_weights.pos')$WGT,PANEL=i)
    
  } else {
    pos_temp<-data.frame(WGT=fread(paste0('/mnt/lustre/groups/biomarkers-brc-mh/TWAS_resource/FUSION/SNP-weights/',i,'/',i,'.pos'))$WGT,PANEL=i)
  }
  pos_all<-rbind(pos_all, pos_temp)
}

# Convert WGT column to match the column names in the predicted expression file
pos_all$ID<-gsub('.*/','',pos_all$WGT)
pos_all$ID<-gsub('.wgt.RDat','',pos_all$ID)
pos_all$ID<-gsub(':','.',pos_all$ID)
pos_all$ID<-gsub('-','.',pos_all$ID)

# Extract features that are to be included in the TWAS
GeneX_incl<-GeneX_all[,(names(GeneX_all) %in% c('FID','IID',pos_all$ID)), with=F]

rm(GeneX_all)
gc()

# Remove features with zero variance (these won't be in the TWAS either)
zeroVar2 <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x)))
  want <- which(out > 1)
  unlist(want)
}

GeneX_incl<-cbind(GeneX_incl[,1:2], GeneX_incl[,-1:-2][,zeroVar2(GeneX_incl[,-1:-2]),with=F])

# See how many genes are being considered
dim(GeneX_incl)[2]-2

# Choose the number of permutations
Nperm<-opt$nperm

# Create progress bar
pb <- txtProgressBar(min = 0, max = Nperm, style = 3)

# Record start time to see how long the permutation procedure takes
start.time <- Sys.time()
start.time

# Each loop is a permutation which generates a random phenotype, test for an association between the random phenotype and every genes. At the end it saves the minimum p-value from each permutation (i.e. null TWAS).
set.seed(opt$seed)
pVal_min<-foreach(j=1:Nperm, .combine=rbind) %dopar% {
  GeneX_incl$Pheno<-rnorm(dim(GeneX_incl)[1])
  
  pVal<-apply(GeneX_incl[,3:(dim(GeneX_incl)[2]-1)], 2, function(x) coef(summary(RcppEigen::fastLm(y=GeneX_incl$Pheno, X=x)))[1,4])
  
  setTxtProgressBar(pb, j)
  print(Sys.time())
  
  if(!is.na(opt$combinations)){
    res<-NULL
    for(i in 1:length(combinations)){
      combinations_i_vec<-unlist(strsplit(combinations[i],' '))
      combination_i_feat<-pos_all[(pos_all$PANEL %in% groups$PANEL[(groups$GROUP %in% combinations_i_vec)]),]
      if(sum((names(pVal) %in% combination_i_feat$ID)) == 0){
        res<-c(res,NA)
      } else {
        res<-c(res,min(pVal[(names(pVal) %in% combination_i_feat$ID)]))
      }
    }
    tmp<-data.frame(t(res))
    names(tmp)<-paste0('combination_',1:length(combinations))
    tmp
  } else {
    tmp<-data.frame(min(pVal))
    tmp
  }
}

# Record when permutation procedure finished and print the time taken.
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Save the minimum p value from each permutation.
write.table(pVal_min, paste0(opt$output,'.Min_P.',Nperm,'_perm.txt'), col.names=T, row.names=F, quote=F)
