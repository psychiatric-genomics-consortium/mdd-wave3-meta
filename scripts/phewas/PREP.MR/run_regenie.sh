# Xshen 30/8/2021


# regenie Step 1 ==============================================================
while read -r a b c d || [ -n "$a" ]; do

jobname="$(echo $(basename ${a}) | sed "s/pheno_//g" | sed "s/.tsv//g")_step1"
fname="${d}_pred.ls"

if [ ! -f ${fname} ]
   then
     qsub -N $jobname scripts/phewas/PREP.MR/regenie_jobs/func.regenie_step1_phewas_gwas.sh "${a}" "${b}" "${c}" "${d}"
     echo $jobname
fi

done < data/MR/regenie_step1_ls.tsv


# regenie Step 2 ==============================================================

# Split files by chromosome
# for chr in {1..22}; do \
#	plink2 --bfile /exports/eddie/scratch/xshen33/ukb_imp_v3.qc --chr $chr --make-bed --out #/exports/eddie/scratch/xshen33/ukb_mdd_gwas/regenie/ukb_imp_v3_${chr}; \
#done

while read -r a b c d || [ -n "$a" ]; do
   jobname="$(echo $(basename ${d}) | sed "s/pheno_//g" )"
   prevjob="$(echo $(basename ${d}) | sed "s/pheno_//g"  | sed "s/_step2/_step1/g")"
   echo $jobname
   echo $prevjob

   for chr in {1..22}; do 
         qsub -N "${jobname}_${chr}" -hold_jid "${prevjob}" scripts/phewas/PREP.MR/regenie_jobs/func.regenie_step2_phewas_gwas.sh $a $b $c $d $chr
   #echo "${jobname}_${chr}"
   done
done < data/MR/regenie_step2_ls.tsv