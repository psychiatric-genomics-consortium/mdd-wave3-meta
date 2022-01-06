#!/bin/sh
#$ -N download_sumstats
#$ -cwd
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=24:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0

Rscript scripts/phewas/PREP.MR/PREP.gwas_sumstats_download.R \
data/MR/MR_pheno_ls_UKB.RData \
data/MR/UKBB_GWAS_Imputed_v3_201807.tsv \
data/MR/BIG40_GWAS.csv \
ata/phe_count.rds \
data/MR/MR_sumstats