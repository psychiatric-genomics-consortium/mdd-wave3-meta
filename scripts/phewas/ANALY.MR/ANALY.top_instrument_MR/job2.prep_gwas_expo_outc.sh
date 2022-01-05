#!/bin/sh
#$ -cwd
#$ -N prep_pheno
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0

Rscript scripts/phewas/ANALY.MR/PREP.pheno_exposure_outcome.R --gwas data/MR/pheno_gwas_forMR.rds \
--ExposureDat data/MR/MR_InterFiles_topSNP/mdd_exposure_outcome/MDD.top_exposure_dat \
--interexp data/MR/MR_InterFiles_topSNP/pheno_exposure_outcome
