#!/bin/sh
#$ -cwd
#$ -N anal_MR
#$ -m beas
#$ -l h_vmem=8G
#$ -l h_rt=6:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0

Rscript scripts/phewas/ANALY.MR/ANALY.top_instrument_MR/ANALY.MR.R --ListMR data/MR/ls.mr_analysis_topSNP.rds \
--out results/phewas/MR/topSNP/ukb/
