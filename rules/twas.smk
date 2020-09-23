##########
# Run twas analysis
##########

####
# Install required software not within conda
####

# install fusion
rule install_fusion:
  output:
    directory("resources/twas/fusion/")
  conda:
    "../envs/twas.yaml"
  shell:
    "git clone git@github.com:gusevlab/fusion_twas.git {output}"

# download plink2R
rule download_plink2R:
  output:
    "resources/twas/plink2R/plink2R-master/data.bed"
  conda:
    "../envs/twas.yaml"
  shell:
    "wget -O resources/twas/plink2R/master.zip https://github.com/gabraham/plink2R/archive/master.zip; unzip resources/twas/plink2R/master.zip -d resources/twas/plink2R"

# install plink2R
rule install_plink2R:
  output:
    directory(".snakemake/conda/6960bedf/lib/R/library/plink2R/")
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript -e 'install.packages(\"resources/twas/plink2R/plink2R-master/plink2R/\",repos=NULL, lib=\".snakemake/conda/6960bedf/lib/R/library/\")'"

# install focus
rule install_focus:
  conda:
    "../envs/twas.yaml"
  output: touch("resources/twas/pyfocus")
  shell:
    "pip install pyfocus==0.6.10 --user"

###
# Munge sumstats
###

# Format sumstats file so FOCUS can read it
rule pre_munge:
  input:
    "results/distribution/daner_pgc_mdd_full_eur_hg19_v3.29.08.gz"
  output:
    "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_premunged.gz"
  conda:
    "../envs/twas.yaml"
  shell:
    "zcat {input} | cut -f 1-5,8-11,19 | sed -e 's/Neff_half/N/g' | gzip -c > {output}"

# munge sumstats using FOCUS munge function
rule focus_munge:
  input:
    premunged="results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_premunged.gz", focus="resources/twas/pyfocus"
  output:
    "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_munged.sumstats.gz"
  conda:
    "../envs/twas.yaml"
  shell:
    "focus munge {input.premunged} --output results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_munged"

###
# Run TWAS
###

# Calculate median effective sample size
rule retrieve_Neff:
  input:
    "results/distribution/daner_pgc_mdd_full_eur_hg19_v3.29.08.gz"
  output:
    "results/twas/median_Neff.txt"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/median_neff.R --daner {input} --out {output}"

# Read in the median Neff
# Neff_file = open("results/twas/median_Neff.txt", "r")
# Neff_char=Neff_file.read()
# Neff_num=float(Neff_char)

# Create list of FUSION SNP-weight sets to be used in the TWAS
weights=["Adrenal_Gland","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Substantia_nigra","CMC.BRAIN.RNASEQ","CMC.BRAIN.RNASEQ_SPLICING","NTR.BLOOD.RNAARR","Pituitary","Thyroid","Whole_Blood","YFS.BLOOD.RNAARR"]

# Create list of chromosome numbers
chr=range(1, 22)

# run twas
rule run_twas:
  input:
    "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_munged.sumstats.gz"
  output:
    "results/twas/PGC_MDD3_twas_{weight}_chr{chr}"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript resources/twas/fusion/FUSION.assoc_test.R "
    "--sumstats {input} "
    "--weights /mnt/lustre/groups/biomarkers-brc-mh/TWAS_resource/FUSION/SNP-weights/{wildcards.weight}/{wildcards.weight}.pos "
    "--weights_dir /mnt/lustre/groups/biomarkers-brc-mh/TWAS_resource/FUSION/SNP-weights/{wildcards.weight} "
    "--ref_ld_chr /scratch/groups/biomarkers-brc-mh/Reference_data/1KG_Phase3/PLINK/EUR/EUR_phase3.MAF_001.chr "
    "--out {output} "
    "--chr {wildcards.chr} "
    "--coloc_P 1e-3 "
    "--GWASN {Neff_num}"

rule fusion_twas:
    input: expand("results/twas/PGC_MDD3_twas_{weight}_chr{chr}", weight=weights, chr=chr)
