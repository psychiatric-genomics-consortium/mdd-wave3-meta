module load igmm/apps/R/3.6.1

rm XS.log
rm -r /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/*

# Prepare exposure data  ----------------------------------------------------------------------

Rscript PREP.exposure_multivar.R --SNPlist ls.snp.multivar_MR.txt \
--CpGlist ls.cpg.multivar_MR.txt \
--mqtl /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/mQTL/mQTL_forMR_GS/ \
--interexp /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/exposure_stats/



# Run MR   ------------------------------------------------------------------------------------
    targetfile="/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/summstats/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta.forPRSice.gz"
    result_tag=`echo $targetfile | sed "s/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta.forPRSice/MDD/"`
    result_tag=`echo $result_tag | sed "s/\\/exports\\/igmm\\/eddie\\/GenScotDepression\\/shen\\/bakup.dat\\/summstats\\///"`  
    result_tag=`echo $result_tag | sed "s/.gz//"`   
    
    rm -r /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/outcome_stats/

    # Prepare outcome data 
 Rscript scripts/phewas/ANALY.MR/ANALY.multi_var/PREP.MDD_outcome.R --exposureSNP data/MR/MR_InterFiles_mvmr/multivar.snps.exposure \
--summout results/distribution/daner_pgc_mdd_noUKBB_eur_hg19_v3.49.24.09.neff.gz \
--interout data/MR/MR_InterFiles_mvmr/

    # MR analysis 
    Rscript scripts/phewas/ANALY.MR/ANALY.multi_var/ANALY.multivar_MR.R \
    --MRexposure data/MR/MR_InterFiles_mvmr/multivar.exposure_dat \
    --MRoutcome data/MR/MR_InterFiles_mvmr/multivar.MDD.outcome_dat \
    --outputfig results/phewas/MR/mvmr/figs \
    --outputtable results/phewas/MR/mvmr/res_multivar.pheno_to_mdd.txt \
    --saveHarmonisedData T
