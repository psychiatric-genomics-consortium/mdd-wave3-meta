#!/bin/sh
#$ -N phewas_xs
#$ -cwd
#$ -m beas
#$ -l h_vmem=16G
#$ -pe sharedmem 3
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0

echo $1

newname=$(echo $1 | sed "s/ /_/g")
echo results/phewas/phewas_out_$newname.rds

Rscript scripts/phewas/ANALY.PRS_phewas/ANALY.phewas.R data/dat.imaging_chunk.rds data/PRS_all.rds results/phewas/models.rds "$1" results/phewas/phewas_out_$newname.rds
