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
  input: rules.download_plink2R.output
  output:
    touch("resources/twas/install_plink2R")
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript -e 'install.packages(\"resources/twas/plink2R/plink2R-master/plink2R/\",repos=NULL)'"

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

# Create list of FUSION SNP-weight sets to be used in the TWAS
weights=["Adrenal_Gland","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Substantia_nigra","CMC.BRAIN.RNASEQ","CMC.BRAIN.RNASEQ_SPLICING","NTR.BLOOD.RNAARR","Pituitary","Thyroid","Whole_Blood","YFS.BLOOD.RNAARR"]

# Create list of chromosome numbers
chr=range(1, 22)

# run twas
rule run_twas:
  input:
    sumstats="results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_munged.sumstats.gz", neff_txt="results/twas/median_Neff.txt", fusion=rules.install_fusion.output, plink2R=rules.install_plink2R.output
  params:
    Neff_num=lambda wildcards, input: float(open(input.neff_txt, "r").read()) 
  output:
    "results/twas/PGC_MDD3_twas_{weight}_chr{chr}"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript resources/twas/fusion/FUSION.assoc_test.R "
    "--sumstats {input.sumstats} "
    "--weights /mnt/lustre/groups/biomarkers-brc-mh/TWAS_resource/FUSION/SNP-weights/{wildcards.weight}/{wildcards.weight}.pos "
    "--weights_dir /mnt/lustre/groups/biomarkers-brc-mh/TWAS_resource/FUSION/SNP-weights/{wildcards.weight} "
    "--ref_ld_chr /scratch/groups/biomarkers-brc-mh/Reference_data/1KG_Phase3/PLINK/EUR/EUR_phase3.MAF_001.chr "
    "--out {output} "
    "--chr {wildcards.chr} "
    "--coloc_P 1e-3 "
    "--GWASN {params.Neff_num}"

rule fusion_twas:
    input: expand("results/twas/PGC_MDD3_twas_{weight}_chr{chr}", weight=weights, chr=chr)

# Run twas using psychENCODE SNP-weights
rule run_psychencode_twas:
  resources: mem_mb=20000 
  input:
    sumstats="results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_munged.sumstats.gz",
    neff="results/twas/median_Neff.txt"
  output:
    "results/twas/psychencode/PGC_MDD3_twas_psychencode_chr{chr}"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript resources/twas/fusion/FUSION.assoc_test.R "
    "--sumstats {input.sumstats} "
    "--weights /scratch/groups/biomarkers-brc-mh/TWAS_resource/PsychEncode/PEC_TWAS_weights/PEC_TWAS_weights.pos "
    "--weights_dir /scratch/groups/biomarkers-brc-mh/TWAS_resource/PsychEncode/PEC_TWAS_weights "
    "--ref_ld_chr /scratch/groups/biomarkers-brc-mh/Reference_data/1KG_Phase3/PLINK/EUR/EUR_phase3.MAF_001.chr "
    "--out {output} "
    "--chr {wildcards.chr} "
    "--coloc_P 1e-3 "
    "--GWASN {Neff_num}"

rule psychencode_twas:
    input: expand("results/twas/psychencode/PGC_MDD3_twas_psychencode_chr{chr}", chr=chrs)

###
# Combine TWAS results
###

rule combine_twas_res:
  input:
    expand("results/twas/PGC_MDD3_twas_{weight}_chr{chr}", weight=weights, chr=chrs), expand("results/twas/psychencode/PGC_MDD3_twas_psychencode_chr{chr}", chr=chrs)
  output:
    "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_TWSig.txt"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/combine_twas.R --out results/twas/twas_results/"

###
# Create Manhattan plots
###

# Manhattan plot based on permutation significance
rule make_manhattan:
  input:
    "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt"
  output:
    "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_Manhattan.png"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript /scratch/groups/biomarkers-brc-mh/TWAS_resource/FUSION/Scripts/Git/opain/TWAS-plotter/TWAS-plotter.V1.0.r \
    --twas {input} \
    --sig_p 1.368572e-06 \
    --output results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_Manhattan \
    --width 5000 \
    --height 5000"

###
# Conditional analysis
###

# Download glist file
rule download_glist:
  output:
    "resources/twas/glist-hg19"
  conda: 
    "../envs/twas.yaml"
  shell:
    "wget -P resources/twas/ https://www.cog-genomics.org/static/bin/plink/glist-hg19"

import pandas as pd
sig_file=pd.read_table("results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_TWSig.txt",sep=' ')
chrs_sig=sig_file.CHR.unique()

# Run conditional analysis
rule run_conditional:
  resources: mem_mb=50000 
  input:
    sig= "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_TWSig.txt",
    ss= "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_munged.sumstats.gz"
  output:
    "results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{chr}.joint_dropped.dat"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript /scratch/groups/biomarkers-brc-mh/TWAS_resource/FUSION/fusion_twas-master/FUSION.post_process.R \
      --input {input.sig} \
      --sumstats {input.ss} \
      --report \
      --ref_ld_chr /scratch/groups/biomarkers-brc-mh/Reference_data/1KG_Phase3/PLINK/EUR/EUR_phase3.MAF_001.chr \
      --out results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{wildcards.chr} \
      --chr {wildcards.chr} \
      --plot \
      --plot_legend all \
      --save_loci \
      --ldsc F \
      --locus_win 500000"

rule conditional:
    input: expand("results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{chr}.joint_dropped.dat", chr=chrs_sig)

######
# Process TWAS results
######
rule process_twas:
  input:
    "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt"
  output:
    "results/twas/twas_results/PGC_MDD3_TWAS_colocalisation.txt"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/process_twas.R"

######
# Process colocalisation analysis results
######
rule process_coloc:
  input:
	"results/twas/twas_results/PGC_MDD3_twas_AllTissues_TWSig_CLEAN.txt"
  output:
    "results/twas/twas_results/PGC_MDD3_TWAS_colocalisation.csv"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/process_coloc.R"


######
# Process conditional analysis results
######

rule process_conditional:
  input: expand("results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{chr}.joint_dropped.dat", chr=chrs_sig), "results/twas/twas_results/PGC_MDD3_twas_AllTissues_TWSig_CLEAN.txt"
  output:
    "results/twas/conditional/PGC_MDD3_TWAS_Conditional_table_novelty.csv"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/process_conditional.R"


########
# Finemap TWAS associations
########

rule run_focus:
  input: "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_munged.sumstats.gz"
  output: "results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv"
  conda: "../envs/twas.yaml"
  shell:
    "focus finemap {input} \
	/scratch/groups/biomarkers-brc-mh/Reference_data/1KG_Phase3/PLINK/EUR/EUR_phase3.MAF_001.chr{wildcards.chr} \
	/scratch/groups/biomarkers-brc-mh/TWAS_resource/FOCUS/MDD_TWAS_db/MDD_TWAS.db \
	--chr {wildcards.chr} \
	--p-threshold 1e-4 \
	--plot \
	--out results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{wildcards.chr}"

rule focus:
    input: expand("results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv", chr=chrs)

rule run_focus_gw:
  input: "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_munged.sumstats.gz"
  output: "results/twas/focus/GW_sig/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv"
  conda: "../envs/twas.yaml"
  shell:
    "focus finemap {input} \
	/scratch/groups/biomarkers-brc-mh/Reference_data/1KG_Phase3/PLINK/EUR/EUR_phase3.MAF_001.chr{wildcards.chr} \
	/scratch/groups/biomarkers-brc-mh/TWAS_resource/FOCUS/MDD_TWAS_db/MDD_TWAS.db \
	--chr {wildcards.chr} \
	--plot \
	--out results/twas/focus/GW_sig/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{wildcards.chr}"

rule focus_gw:
    input: expand("results/twas/focus/GW_sig/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv", chr=chrs)

# Process the FOCUS results
rule process_focus:
  input: expand("results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv", chr=chrs)
  output:
    "results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.results.csv"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/process_focus.R"

########
# TWAS-GSEA
########

# Split twas results by tissue type
rule prep_for_gsea:
  input: "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt"
  output:
    "results/twas/twas_results/PGC_MDD3_twas_BLOOD_GW.txt"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/prep_for_gsea.R"

# Run twas-gsea across results
rule run_twas_gsea:
  resources: mem_mb=30000, cpus=3
  input: ss="results/twas/twas_results/PGC_MDD3_twas_{set}_GW.txt", gmt="resources/twas/gsea/c2_and_c5.all.v7.1.entrez.gmt"
  output:
    "results/twas/gsea/PGC_MDD3_twas_{set}_GSEA.log"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript /mnt/lustre/groups/biomarkers-brc-mh/TWAS_resource/FUSION/Scripts/Git/opain/TWAS-GSEA/TWAS-GSEA.V1.2.R \
  --twas_results {input.ss} \
  --pos results/twas/twas_results/PGC_MDD3_twas.pos \
  --gmt_file {input.gmt} \
  --qqplot F \
  --expression_ref /users/k1806347/brc_scratch/Analyses/Lorenza/Clean/TWAS-GSEA/FUSION_PsychENCODE_FeaturePredictions.csv.gz \
  --self_contained F \
  --min_r2 0.05 \
  --linear_p_thresh 1 \
  --n_cores 3 \
  --competitive T \
  --covar GeneLength,NSNP \
  --output results/twas/gsea/PGC_MDD3_twas_{wildcards.set}_GSEA"

rule twas_gsea:
    input: expand("results/twas/gsea/PGC_MDD3_twas_{set}_GSEA.log", set=["AllTissues","BRAIN","BLOOD","HPA","HPT"])

########
# Run MESC
########

# Estimate h2med
rule est_h2med:
  input: "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_munged.sumstats.gz"
  output: "results/twas/mesc/PGC_MDD3_TWAS.MESC.all.h2med"
  conda: "../envs/mesc.yaml"
  shell: "resources/twas/mesc/mesc/run_mesc.py --h2med {input} --exp-chr resources/twas/mesc/GTEx_v8/All_Tissues --out results/twas/mesc/PGC_MDD3_TWAS.MESC"

# MESC estimates non-significant SNP-based heritability mediated via cis-regulated expression. This is a surprising given the many colocalised signals identified in the TWAS, and it makes question the validity of estimates provided by MESC.

# Estimate gene-set enrichment
rule est_h2med_sets:
  resources: mem_mb=20000 
  input: "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.29.08_munged.sumstats.gz"
  output: "results/twas/mesc/PGC_MDD3_TWAS.MESC.sets.all.h2med"
  conda: "../envs/mesc.yaml"
  shell: "resources/twas/mesc/mesc/run_mesc.py --h2med {input} --exp-chr resources/twas/mesc/GTEx_gene_set_expscores/All_Tissues --out results/twas/mesc/PGC_MDD3_TWAS.MESC.sets"
