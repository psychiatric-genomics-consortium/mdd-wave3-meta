IMPV3=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3
SCRATCH=/exports/eddie/scratch/$USER

# European ancestries analysis

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--pfile $IMPV3/pgen/ukb_imp_v3.autosome.qc.mac100.info9.hwe10 \
--keep /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/v2/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id \
--extract $1 \
--maf 0.01 \
--mac 100 \
--geno 0.02 \
--mind 0.02 \
--make-bed \
--out $2 \
--memory 15360 \
--threads 4