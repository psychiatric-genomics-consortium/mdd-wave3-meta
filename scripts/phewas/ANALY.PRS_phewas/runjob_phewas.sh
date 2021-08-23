while read targetfile; do

fname="results/phewas/phewas_out_$(echo $targetfile | sed "s/ /_/g").rds"
jobname=$(echo $targetfile | sed "s/ /_/g")

if [ ! -f ${fname} ]
   then
     qsub scripts/phewas/ANALY.PRS_phewas/job.phewas.sh "${targetfile}" 
     echo $fname
fi

done < data/phewas_categories.tsv