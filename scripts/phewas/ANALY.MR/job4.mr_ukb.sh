#!/bin/sh
#$ -cwd
#$ -N prep_MR
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=6:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/3.6.1

Rscript scripts/phewas/ANALY.MR/ANALY.MR.R --ListMR data/MR/ls.mr_analysis.rds \
--out results/phewas/MR/ukb/
