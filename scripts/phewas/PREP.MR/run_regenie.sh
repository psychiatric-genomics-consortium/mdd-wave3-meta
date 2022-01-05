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

while read -r a b c d || [ -n "$a" ]; do
   jobname="$(echo $(basename ${d}) | sed "s/pheno_//g" )"
   prevjob="$(echo $(basename ${d}) | sed "s/pheno_//g"  | sed "s/_step2/_step1/g")"
   echo $jobname
   echo $prevjob

   qsub -N "${jobname}" -hold_jid "${prevjob}" scripts/phewas/PREP.MR/regenie_jobs/func.regenie_step2_phewas_gwas.sh $a $b $c $d

done < data/MR/regenie_step2_ls.tsv


# regenie Step 2 (remaining) ==================================================

while read -r a b c d f || [ -n "$a" ]; do
   jobname="$(echo $(basename ${d}) | sed "s/pheno_//g" )"
   echo $jobname
   echo "starting block:" $f

   qsub -N "${jobname}" scripts/phewas/PREP.MR/regenie_jobs/func.regenie_step2_phewas_gwas_remaining.sh $a $b $c $d $f

done < data/MR/regenie_step2_ls_remaining.tsv
