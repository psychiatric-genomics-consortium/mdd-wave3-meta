#!/bin/sh
#$ -cwd
#$ -N prep_mdd_outcome
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 2
#$ -l h_rt=6:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0

Rscript scripts/phewas/ANALY.MR/ANALY.Howard_etal_MR/PREP.MDD_outcome.R --exposure data/MR/MR_InterFiles_Howard/pheno_exposure_outcome \
--summout results/distribution/Howard_mdd_PGC23andme.gz \
--interout data/MR/MR_InterFiles_Howard/mdd_exposure_outcome/