#!/bin/sh
#$ -N reformat_gwas
#$ -cwd
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=5:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0


Rscript scripts/phewas/PREP.MR/PREP.ukb_gwas_sumstats_reformat.R /exports/eddie/scratch/xshen33/phewas_gwas data/MR/MR_sumstats