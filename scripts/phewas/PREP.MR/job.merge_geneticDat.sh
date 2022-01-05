#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -cwd


cd /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/

# rm -f list_pgen.txt
# ls /exports/eddie/scratch/xshen33/ukb_mdd_gwas/regenie/*.pgen | sed "s/.pgen//g" > list_pgen.txt
# ls /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/*.pgen | sed "s/\\.pgen//g" >> list_pgen.txt

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--pmerge-list list_pgen.txt \
--make-pgen \
--out ukb_imp_v3.qc \
--memory 22000 \
--threads 12

