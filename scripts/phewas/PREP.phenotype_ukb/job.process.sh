#!/bin/sh
#$ -N phenotype_correction
#$ -cwd
#$ -m beas
#$ -l h_vmem=64G
#$ -l h_rt=3:00:00
. /etc/profile.d/modules.sh

module load igmm/apps/R/3.6.1

Rscript scripts/phewas/PREP.phenotype_ukb/process_BrainImaging_2.R data/Data_Dictionary_Showcase.csv results/phewas/data_dictionary/fields.imaging_phenotype.txt data/2021-04-phenotypes-ukb44797/Imaging.rds results/phewas/data_dictionary/fields.brain_imaging_QC_cov.txt data/2021-04-phenotypes-ukb44797/Recruitment.rds data/dat.imaging_chunk.rds results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt

Rscript scripts/phewas/PREP.phenotype_ukb/process_loose_fields_7.R data/loose_fields data/Codings.csv data/2021-04-phenotypes-ukb44797/Imaging.rds results/phewas/data_dictionary/category_loose_fields data/dat.loose_field_chunk.rds results/phewas/data_dictionary/fields.final.loose_field.txt