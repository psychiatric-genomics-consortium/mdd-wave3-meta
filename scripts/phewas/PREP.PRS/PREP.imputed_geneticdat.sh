#$ -N target_qc
#$ -l h_vmem=4G
# -pe sharedmem 6
#$ -l h_rt=4:00:00
#$ -cwd

# Example job to QC UKB imputed PGEN data for PRS target 

IMPV3=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3
SCRATCH=/exports/eddie/scratch/$USER

# European ancestries analysis

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--pfile $IMPV3/pgen/ukb_imp_v3.autosome.qc.mac100.info9.hwe10 \
--keep /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/v2/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id \
--extract ../../../data/mdd3.snps \
--maf 0.005 \
--mac 100 \
--geno 0.02 \
--mind 0.02 \
--make-bed \
--out $SCRATCH/ukb_imp_v3.qc \
--memory 15360 \
--threads 6