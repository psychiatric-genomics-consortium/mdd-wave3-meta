#!/bin/sh
#$ -N pgc3_mddprs
#$ -cwd
#$ -m beas
#$ -l h_vmem=64G
#$ -pe sharedmem 5
#$ -l h_rt=12:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/3.6.1

Rscript /exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice.R \
    --prsice /exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice_linux \
    --base /exports/eddie/scratch/xshen33/23andmePGC_June13_GCldsc1.meta_forPRSice \
    --target /exports/eddie/scratch/xshen33/ukb_imp_v3.qc \
    --thread 5 \
    --stat BETA \
    --pheno data/ukb_dummypheno \
    --pheno-col dummy_pheno \
    --clump-kb 250 \
    --clump-r2 0.25 \
    --bar-levels 0.00000005,0.0000001,0.000001,0.00001,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
    --fastscore \
    --all-score \
    --out data/mdd_prs_imputed
