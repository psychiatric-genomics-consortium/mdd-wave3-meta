#!/bin/sh
#$ -cwd
#$ -N prep_MR
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/3.6.1

Rscript scripts/phewas/ANALY.MR/PREP.pheno_exposure_outcome.R --gwas data/MR/pheno_gwas_forMR.rds \
--ExposureDat data/MR/MR_InterFiles/mdd_exposure_outcome/MDD.exposure_dat \
--interexp data/MR/MR_InterFiles/pheno_exposure_outcome