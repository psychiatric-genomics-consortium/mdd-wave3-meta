#!/bin/sh
#$ -N pgc3_mddprs
#$ -cwd
#$ -m beas
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/3.6.1

Rscript /exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice.R \
    --prsice /exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice_linux \
    --base ../../data/mdd_noUKB_forPRSice \
    --target /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_genotype/autosome.qc.maf01.hwe5e-6.geno02.mind02.snps \
    --thread 5 \
    --stat OR \
    --pheno ../../data/ukb_dummypheno \
    --pheno-col dummy_pheno \
    --clump-kb 250 \
    --clump-r2 0.25 \
    --bar-levels 0.00000005,0.0000001,0.000001,0.00001,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
    --fastscore \
    --all-score \
    --out ../../data/mdd_prs_genotyped
