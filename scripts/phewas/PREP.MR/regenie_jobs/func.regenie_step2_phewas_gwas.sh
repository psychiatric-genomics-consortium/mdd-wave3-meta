#!/bin/sh
#$ -cwd
#$ -m beas
#$ -l h_vmem=4G
#$ -pe sharedmem 4
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh


/exports/igmm/eddie/GenScotDepression/shen/Tools/regenie/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
  --step 2 \
  --bed /exports/eddie/scratch/xshen33/ukb_mdd_gwas/regenie/ukb_imp_v3_${5} \
  --ref-first \
  --phenoFile ${1} \
  --covarFile ${2} \
  --firth 0.01 \
     --approx \
  --pred ${3} \
  --bsize 400 \
  --threads 10 \
    --no-split \
     --gz \
  --out ${4}_${5}
