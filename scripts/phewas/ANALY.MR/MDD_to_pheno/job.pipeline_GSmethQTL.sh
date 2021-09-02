
rm XS.log
rm -r data/MR/MR_InterFiles/*

# Prepare exposure data  ----------------------------------------------------------------------

Rscript scripts/phewas/ANALY.MR/PREP.MDD_exposure.R --exposure results/distribution/daner_pgc_mdd_noUKBB_eur_hg19_v3.49.24.05.rp.gz \
--interexp data/MR/MR_InterFiles/exposure_stats/



# Run MR   ------------------------------------------------------------------------------------

while read targetfile; do
    result_tag=`echo $targetfile | sed "s/.rds//"`
    result_tag=`echo $result_tag | sed "s/\\/exports\\/igmm\\/eddie\\/GenScotDepression\\/shen\\/ActiveProject\\/Genetic\\/MR_meth_MDD\\/data\\/mQTL\\/mQTL_forMR_GS\\/mQTL.//"`        
    
    rm -r /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/MR_InterFiles/outcome_stats/

    # Prepare outcome data 
    Rscript PREP.outcome.R --exposure /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/MR_InterFiles/exposure_stats/ \
    --summout ${targetfile} \
    --outcomeName ${result_tag} \
    --interout /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/MR_InterFiles/outcome_stats/

    # MR analysis 
    Rscript ANALY.MR.R \
    --MRexposure /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/MR_InterFiles/exposure_stats/ \
    --MRoutcome /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/MR_InterFiles/outcome_stats/ \
    --outcomeName ${result_tag} \
    --outputfig /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/result/GS_MR_MDD_to_DNAm/MDD_to_${result_tag}_figs/ \
    --outputtable /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/result/GS_MR_MDD_to_DNAm/Summary_MDD_to_${result_tag}.csv \
    --saveHarmonisedData T
done < ls.outcome

