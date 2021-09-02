#!/bin/sh
#$ -cwd
#$ -m beas
#$ -l h_vmem=16G
#$ -pe sharedmem 4
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh


/exports/igmm/eddie/GenScotDepression/shen/Tools/regenie/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
  --step 1 \
  --bed /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_genotype/autosome.qc.maf01.hwe5e-6.geno02.mind02.snps \
  --phenoFile $1 \
  --covarFile $2 \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix $3 \
  --out $4