#!/bin/sh
#$ -N reformat_gwas
#$ -cwd
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=10:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0

Rscript scripts/phewas/PREP.MR/PREP.gwas_sumstats_format_for_meta.R \
data/MR/MR_pheno_ls_UKB.RData \
data/MR/meta.ls.rds \
data/MR/MR_sumstats/meta \
data/MR/UKBB_GWAS_Imputed_v3_201807.tsv \
data/MR/BIG40_GWAS.csv \
data/MR/MR_sumstats