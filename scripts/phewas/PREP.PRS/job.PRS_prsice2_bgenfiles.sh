#!/bin/sh
#$ -N pgc3_mddprs
#$ -cwd
#$ -m beas
#$ -l h_vmem=32G
#$ -pe sharedmem 8
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/3.5.1

Rscript /exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice.R \
    --prsice /exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice_linux \
    --base /exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/data/daner_pgc_mdd_noUKBB_eur_hg19_v3.229.09_forPRS.txt \
    --target '/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3/ukb_imp_chr#_v3' \
    --type bgen \
    --maf 0.005 \
    --keep /exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/data/ukb_unrelated_white_ID.txt \
    --extract /exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/data/mdd_prs/bgen_PRS/pgc3_meta_ukb_PRS.valid \
    --thread 30 \
    --stat OR \
    --binary-target T,T,T,T \
    --pheno /exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/data/ukb_mdd_phenotypes_bgen.txt \
    --pheno-col mdd_nerves,mdd_smith,mdd_icd,mdd_cidi \
    --dose-thres 0.9 \
    --clump-kb 250 \
    --clump-r2 0.25 \
    --bar-levels 0.00000005,0.0000001,0.000001,0.00001,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
    --fastscore \
    --all-score \
    --allow-inter \
    --out /exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/data/mdd_prs/bgen_PRS/pgc3_meta_ukb_PRS
