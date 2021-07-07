#!/bin/sh
#$ -N MEGApgrs
#$ -cwd
#$ -m beas
#$ -l h_vmem=128G
#$ -l h_rt=24:00:00
. /etc/profile.d/modules.sh

module add igmm/apps/R/3.2.4


R --file=/exports/igmm/eddie/GenScotDepression/shen/PGRS/PRSice_v1.25/PRSice_v1.25.R -q --args \
plink /exports/igmm/eddie/GenScotDepression/shen/PGRS/plink_linux_x86_64/plink \
base /exports/igmm/eddie/GenScotDepression/shen/PGRS/mega_June132018_MDD_PGR/MDD_meta_23andmePGCnoMRIUKB_20181029_PGRS/base.dat/ii.basedat.MDD.3cohortsmeta \
target /exports/igmm/eddie/GenScotDepression/shen/PGRS/mega_June132018_MDD_PGR/MDD_meta_23andmePGCnoMRIUKB_20181029_PGRS/target.dat/autosome.qc.maf01.hwe5e-6.geno02.mind02.snps \
clump.p1 1 \
clump.p2 1 \
clump.r2 0.25 \
clump.kb 250 \
fastscore T \
barchart.levels 0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
report.individual.scores T \
report.best.score.only F \
covary F \
no.regression F \
debug.mode T \
plink.silent F \
binary.target T \
figname megaPGRS \
pheno.file /exports/igmm/eddie/GenScotDepression/shen/PGRS/mega_June132018_MDD_PGR/MDD_meta_23andmePGCnoMRIUKB_20181029_PGRS/dummy.pheno.MDDpgrs \
wd /exports/igmm/eddie/GenScotDepression/shen/PGRS/mega_June132018_MDD_PGR/MDD_meta_23andmePGCnoMRIUKB_20181029_PGRS/MDD_PGRS_23andmePGCnoMRIUKB_20181029.XS/prsice.output/

