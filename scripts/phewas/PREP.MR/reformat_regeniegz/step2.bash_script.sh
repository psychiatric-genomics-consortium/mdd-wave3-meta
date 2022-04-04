#!/bin/sh
#$ -cwd
#$ -m beas
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0


# This script is used to parse gz format sumstats from regenie by each phenotype

cd /exports/eddie/scratch/xshen33/phewas_gwas/
bash reformat_singlerun.sh