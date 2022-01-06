#!/bin/sh
#$ -cwd
#$ -N prep_MR
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=6:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/3.6.1

Rscript scripts/phewas/ANALY.MR/ANALY.MR_IEU.R --ListMR data/MR/ls.mr_analysis.rds \
--MDDsumstats results/distribution/daner_pgc_mdd_noUKBB_eur_hg19_v3.49.24.05.rp.gz \
--MDDinstrument data/MR/MR_InterFiles/mdd_exposure_outcome/MDD.exposure_dat \
--out results/phewas/MR/mr_base/
