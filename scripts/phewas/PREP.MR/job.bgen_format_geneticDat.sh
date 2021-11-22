#$ -l h_rt=48:00:00
#$ -l h_vmem=2G
#$ -pe sharedmem 12
#$ -t 16-22
#$ -tc 8
#$ -e logs
#$ -o logs
#$ -cwd

CHR=$SGE_TASK_ID

IMPV3=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3
SCRATCH=/exports/eddie/scratch/$USER/ukb_mdd_gwas/regenie

mkdir -p $SCRATCH

cat $IMPV3/ukb_mfi_chr${CHR}_v3.txt | awk '{if($6 >= 0.005 && $8 >= 0.1) {print $2}}' > $SCRATCH/ukb_chr${CHR}.qc.snps
wc -l $SCRATCH/ukb_chr${CHR}.qc.snps

/exports/igmm/eddie/GenScotDepression/local/bin/bgenix1.1.4 \
-g $IMPV3/ukb_imp_chr${CHR}_v3.bgen \
-incl-rsids $SCRATCH/ukb_chr${CHR}.qc.snps > $SCRATCH/ukb_imp_chr${CHR}_v3.qc.bgen

cp $IMPV3/ukb4844_imp_chr*.sample $SCRATCH
