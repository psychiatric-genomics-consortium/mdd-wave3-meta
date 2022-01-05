#$ -l h_rt=48:00:00
#$ -l h_vmem=2G
#$ -pe sharedmem 12
#$ -t 3
#$ -tc 3
#$ -e logs
#$ -o logs
#$ -cwd

CHR=$SGE_TASK_ID

IMPV3=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3
SCRATCH=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen

mkdir -p $SCRATCH

### create snp lists (info>0.1, MAF>0.001)

cat $IMPV3/ukb_mfi_chr${CHR}_v3.txt | awk '{if($6 >= 0.001 && $8 >= 0.1) {print $2}}' > $SCRATCH/ukb_chr${CHR}.qc.snps
wc -l $SCRATCH/ukb_chr${CHR}.qc.snps

/exports/igmm/eddie/GenScotDepression/local/bin/bgenix1.1.4 \
-g $IMPV3/ukb_imp_chr${CHR}_v3.bgen \
-incl-rsids $SCRATCH/ukb_chr${CHR}.qc.snps > $SCRATCH/ukb_imp_chr${CHR}_v3.qc.bgen


### make pgen files per chr (mac>20, no dup, no multi-allelic)

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--bgen $SCRATCH/ukb_imp_chr${CHR}_v3.qc.bgen ref-first \
--sample $IMPV3/ukb4844_imp_chr${CHR}_v3_s487395.sample \
--mac 20 \
--rm-dup 'force-first' \
--max-alleles 2 \
--make-pgen \
--out $SCRATCH/ukb_imp_chr${CHR}_v3.qc \
--memory 22000 \
--threads 12

rm $SCRATCH/ukb_imp_chr${CHR}_v3.qc.bgen


### merge pgen files (job.merge_geneticDat.sh)


