#!/bin/sh
#$ -cwd
#$ -m beas
#$ -l h_vmem=2G
#$ -pe sharedmem 15
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh


/exports/igmm/eddie/GenScotDepression/shen/Tools/regenie/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
  --step 2 \
  --pgen /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/ukb_imp_v3.qc \
  --ref-first \
  --phenoFile ${1} \
  --covarFile ${2} \
  --firth 0.01 \
  --approx \
  --pred ${3} \
  --bsize 400 \
  --threads 32 \
  --no-split \
  --starting-block ${5} \
  --gz \
  --out ${4}
