#!/bin/sh
#$ -cwd
#$ -N prep_mdd_outcome_topsnp
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=6:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0

Rscript scripts/phewas/ANALY.MR/ANALY.top_instrument_MR/PREP.MDD_outcome.R --exposure data/MR/MR_InterFiles_topSNP/pheno_exposure_outcome \
--summout results/distribution/daner_pgc_mdd_noUKBB_eur_hg19_v3.49.24.09.neff.gz \
--interout data/MR/MR_InterFiles_topSNP/mdd_exposure_outcome/