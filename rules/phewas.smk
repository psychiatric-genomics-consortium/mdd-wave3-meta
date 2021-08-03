# Environment: R>3.6, Linux
# R packages: data.table, dplyr, tidyverse, optparse
# Other: PRSice 2.0
# Temporary data stored in the data folder

# PRS  =============================================================================

### PREP for PRSice
# Prep sumstats
rule reformat_sumstats_PRSice:
    input:
        "results/distribution/daner_pgc_mdd_noUKBB_eur_hg19_v3.49.24.05.rp.gz"
    output:
        "data/mdd_noUKB_forPRSice"
    shell:
        "Rscript script/phewas/PREP.PRS/PREP.sumstats_forPRSice.R --sumstats {input} --out {output}"

# Prep dummy phenotype
rule dummypheno_PRSice:
    input:
        "/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_genotype/autosome.qc.maf01.hwe5e-6.geno02.mind02.snps.fam",
                    "/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/v2/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id"
    output:
        "data/ukb_dummypheno"
    shell:
        "Rscript script/phewas/PREP.PRS/PREP.pheno_forPRSice.R {input} {output}"

# Run PRSice on genotyped data (or: qsub script/phewas/PREP.PRS/job.PRS_genotyped.sh)
rule mddprs_genotyped:
    input:
        "data/mdd_noUKB_forPRSice"
        "/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_genotype/autosome.qc.maf01.hwe5e-6.geno02.mind02.snps",
        "data/ukb_dummypheno"
    output:
        "data/mdd_prs_genotyped"
    shell:
        "bash script/phewas/PREP.PRS/job.PRS_genotyped_snkmk.sh {input} {output}"

# Prepare snp list for imputed data 
rule mksnpls_imputed:
    input:
        "data/mdd_noUKB_forPRSice"
    output:
        "data/mdd3.snps"
    shell:
        "cat {input} | tail -n +2 | awk '{{print $2}}' > {output}"

# Prepare imputed genetic data in plink format: Extract SNPs from PGEN file, QC, and convert to BED format (or: qsub script/phewas/PREP.PRS/PREP.imputed_geneticdat.sh)
# Check: https://github.com/ccbs-stradl/ukb_impv_prs
rule bedformat_imputed:
    input:
        "data/mdd3.snps"
    output:
        "\$SCRATCH/ukb_imp_v3.qc"
    shell:
        "bash script/phewas/PREP.PRS/PREP.imputed_geneticdat_snkmk.sh {input} {output}"

# Run PRSice on imputed data (or: qsub script/phewas/PREP.PRS/job.PRS_imputed.sh)
# Check: https://github.com/ccbs-stradl/ukb_impv_prs
rule mddprs_imputed:
    input:
        "data/mdd_noUKB_forPRSice"
        "\$SCRATCH/ukb_imp_v3.qc",
        "data/ukb_dummypheno"
    output:
        "data/mdd_prs_genotyped"
    shell:
        "bash script/phewas/PREP.PRS/job.PRS_imputed_snkmk.sh {input} {output}"


# PheWAS data prep =====================================================================
# Download data dictionary
rule ukb_data_dictionary:
    input:
         HTTP.remote("https://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.csv")
    output:
        "data/Data_Dictionary_Showcase.csv"
    shell:
        "cp {input} {output}"

rule data_coding:
    input:
         HTTP.remote("https://biobank.ctsu.ox.ac.uk/~bbdatan/Codings.csv")
    output:
        "data/Codings.csv"
    shell:
        "cp {input} {output}"

rule loose_pheno:
    input:
        "data/Data_Dictionary_Showcase.csv",
        "results/phewas/data_dictionary/field_experiment_procedure.txt",
        "results/phewas/data_dictionary/field_covered_by_summary_scores.txt"
    output:
        "data/loose_fields"
    shell:
        "Rscript scripts/phewas/PREP.phenotype_ukb/fields_general_1.R {input} {output}"

rule chunk_brain_imaging:
    input:
        "data/Data_Dictionary_Showcase.csv",
        "results/phewas/data_dictionary/fields.imaging_phenotype.txt",
        "data/2021-04-phenotypes-ukb44797/Imaging.rds",
        "results/phewas/data_dictionary/fields.brain_imaging_QC_cov.txt",
        "data/2021-04-phenotypes-ukb44797/Recruitment.rds"
    output:
        "data/dat.imaging_chunk.rds",
        "results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt"
    shell:
        "Rscript scripts/phewas/PREP.phenotype_ukb/process_BrainImaging_2.R {input} {output}"

rule chunk_cognition:
    input:
        "data/2021-04-phenotypes-ukb44797/CognitiveFunction.rds",
        "data/2021-04-phenotypes-ukb44797/CognitiveFunctionOnline.rds",
        "results/phewas/data_dictionary/fields.cog_after_manual_check.txt"
    output:
        "data/dat.cognition_chunk.rds",
        "results/phewas/data_dictionary/fields.final.cognition.txt"
    shell:
        "Rscript scripts/phewas/PREP.phenotype_ukb/process_cognition_3.R {input} {output}"

rule chunk_diet:
    input:
        "data/Data_Dictionary_Showcase.csv",
        "data/Codings.csv",
        "data/2021-04-phenotypes-ukb44797/DietByHourRecall.rds",
        "data/2021-04-phenotypes-ukb44797/Touchscreen.rds"
    output:
        "data/dat.diet_chunk.rds",
        "results/phewas/data_dictionary/fields.final.diet.txt"
    shell:
        "Rscript scripts/phewas/PREP.phenotype_ukb/process_diet_4.R {input} {output}"

rule chunk_activity:
    input:
        "data/Data_Dictionary_Showcase.csv",
        "data/Codings.csv",
        "data/2021-04-phenotypes-ukb44797/PhysicalActivityMeasurement.rds",
    output:
        "data/dat.activity_chunk.rds",
        "results/phewas/data_dictionary/fields.final.activity_data_and_QC.txt"
    shell:
        "Rscript scripts/phewas/PREP.phenotype_ukb/process_PhysicalActivity_5.R {input} {output}"

