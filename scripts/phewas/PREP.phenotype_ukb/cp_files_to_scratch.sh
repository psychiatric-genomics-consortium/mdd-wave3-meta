# log on a staging node 
# qlogin -q staging

# Latest release
cp -r /exports/igmm/datastore/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797 $myscratch
# Chronic pain
cp /exports/igmm/datastore/GenScotDepression/data/ukb/phenotypes/fields/2021-04-chronicpain-ukb45596/* $myscratch/2021-04-phenotypes-ukb44797/
# Blood biomarkers
cp /exports/igmm/datastore/GenScotDepression/data/ukb/phenotypes/fields/2020-02-imaging-ukb40531/* $myscratch/2021-04-phenotypes-ukb44797/

# Other derived data
mkdir $myscratch/2021-04-phenotypes-ukb44797/derived
# MHQ
cp /exports/igmm/datastore/GenScotDepression/data/ukb/phenotypes/psych/mhq/MHQ.1907* $myscratch/2021-04-phenotypes-ukb44797/derived/
# MDD definitions
cp /exports/igmm/datastore/GenScotDepression/data/ukb/phenotypes/mood/mdd/mdd_pipeline/ukb24262-mdd.mdd_phenotypes.rds $myscratch/2021-04-phenotypes-ukb44797/derived/
# Genetic covariates for people eligible for PRS analysis
cp /exports/igmm/datastore/GenScotDepression/data/ukb/genetics/gwas/covar/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.covars.rds $myscratch/2021-04-phenotypes-ukb44797/derived/

# Make a hyperlink
ln -s /exports/eddie/scratch/xshen33/2021-04-phenotypes-ukb44797/ data/2021-04-phenotypes-ukb44797/

