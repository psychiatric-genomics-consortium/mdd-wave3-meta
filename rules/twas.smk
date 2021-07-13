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

# install TWAS-plotter
rule install_twas_plotter:
  output:
    directory("resources/twas/TWAS-plotter/")
  conda:
    "../envs/twas.yaml"
  shell:
    "git clone git@github.com:opain/TWAS-plotter.git {output}"

# install TWAS-GSEA
rule install_twas_gsea:
  output:
    directory("resources/twas/TWAS-GSEA/")
  conda:
    "../envs/twas.yaml"
  shell:
    "git clone git@github.com:opain/TWAS-GSEA.git {output}"

# Download gene-sets
rule download_c2_sets:
  output:
    "resources/twas/gsea/c2.all.v7.4.entrez.gmt"
  conda:
    "../envs/twas.yaml"
  shell:
    "wget -O resources/twas/gsea/c2.all.v7.4.entrez.gmt --no-check-certificate 'https://docs.google.com/uc?export=download&id=1msbXzWTqBNa3pxVqzOwCjb7SSgpXykik'"
    
rule download_c5_sets:
  output:
    "resources/twas/gsea/c5.all.v7.4.entrez.gmt"
  conda:
    "../envs/twas.yaml"
  shell:
    "wget -O resources/twas/gsea/c5.all.v7.4.entrez.gmt --no-check-certificate 'https://docs.google.com/uc?export=download&id=1lzpZP4TtZjIr94rKjNbwWQYLgH7u70gL'"
    
# Combine gene-sets into a single file
rule combine_gene_sets:
  input:
    c2="resources/twas/gsea/c2.all.v7.4.entrez.gmt", c5="resources/twas/gsea/c5.all.v7.4.entrez.gmt"
  output:
    "resources/twas/gsea/c2_and_c5.all.v7.4.entrez.gmt"
  conda:
    "../envs/twas.yaml"
  shell:
    "cat {input.c2} {input.c5} > resources/twas/gsea/c2_and_c5.all.v7.4.entrez.gmt"

# Install MESC
rule install_mesc:
  output:
    directory("resources/twas/mesc/")
  conda:
    "../envs/twas.yaml"
  shell:
    "git clone git@github.com:douglasyao/mesc.git {output}"

# Download MESC expression scores
rule download_mesc_score:
  input:
    mesc_dir=rules.install_mesc.output
  output:
    directory("resources/twas/mesc_data/GTEx_v8/")
  conda:
    "../envs/twas.yaml"
  shell:
    "wget -O resources/twas/mesc_data/All_Tissues.tar.gz http://gusevlab.org/projects/MESC/All_Tissues.tar.gz ; tar -xvf resources/twas/mesc_data/All_Tissues.tar.gz -C resources/twas/mesc_data/; rm resources/twas/mesc_data/All_Tissues.tar.gz"

# Download MESC expression scores for sets
rule download_mesc_score_sets:
  input:
    mesc_dir=rules.install_mesc.output
  output:
    directory("resources/twas/mesc_data/GTEx_gene_set_expscores")
  conda:
    "../envs/twas.yaml"
  shell:
    "wget -O resources/twas/mesc_data/All_Tissues_Gene_Sets.tar.gz http://gusevlab.org/projects/MESC/All_Tissues_Gene_Sets.tar.gz ; tar -xvf resources/twas/mesc_data/All_Tissues_Gene_Sets.tar.gz -C resources/twas/mesc_data/; rm resources/twas/mesc_data/All_Tissues_Gene_Sets.tar.gz"
    
###
# Munge sumstats
###

# Format sumstats file so FOCUS can read it
rule pre_munge:
  input:
    "results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp.gz"
  output:
    "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_premunged.gz"
  conda:
    "../envs/twas.yaml"
  shell:
    "zcat {input} | cut -f 1-5,8-11,19 | sed -e 's/Neff_half/N/g' | gzip -c > {output}"

# munge sumstats using FOCUS munge function
rule focus_munge:
  input:
    premunged="results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_premunged.gz", focus="resources/twas/pyfocus"
  output:
    "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_munged.sumstats.gz"
  conda:
    "../envs/twas.yaml"
  shell:
    "focus munge {input.premunged} --output results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_munged"

###
# Run TWAS
###

# Calculate median effective sample size
rule retrieve_Neff:
  input:
    "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_munged.sumstats.gz"
  output:
    "results/twas/median_Neff.txt"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/median_neff.R --munged {input} --out {output}"

# Create list of FUSION SNP-weight sets to be used in the TWAS
weights=["Adrenal_Gland","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Substantia_nigra","CMC.BRAIN.RNASEQ","CMC.BRAIN.RNASEQ_SPLICING","NTR.BLOOD.RNAARR","Pituitary","Thyroid","Whole_Blood","YFS.BLOOD.RNAARR"]

# Create list of chromosome numbers
chr=range(1, 22)

# run twas
rule run_twas:
  resources: mem_mb=20000
  input:
    sumstats="results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_munged.sumstats.gz", neff_txt="results/twas/median_Neff.txt", fusion=rules.install_fusion.output, plink2R=rules.install_plink2R.output
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
    sumstats="results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_munged.sumstats.gz", neff_txt="results/twas/median_Neff.txt", fusion=rules.install_fusion.output, plink2R=rules.install_plink2R.output
  params:
    Neff_num=lambda wildcards, input: float(open(input.neff_txt, "r").read()) 
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
    "--GWASN {params.Neff_num}"

rule psychencode_twas:
    input: expand("results/twas/psychencode/PGC_MDD3_twas_psychencode_chr{chr}", chr=chr)

###
# Combine TWAS results
###

rule combine_twas_res:
  input:
    expand("results/twas/PGC_MDD3_twas_{weight}_chr{chr}", weight=weights, chr=chr), expand("results/twas/psychencode/PGC_MDD3_twas_psychencode_chr{chr}", chr=chr)
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
    results="results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt", fusion=rules.install_twas_plotter.output
  output:
    "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_Manhattan.png"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript resources/twas/TWAS-plotter/TWAS-plotter.V1.0.r \
    --twas {input.results} \
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

# Run conditional analysis
sig_file = Path("results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_TWSig.txt")
if sig_file.is_file():
  import pandas as pd
  sig_df=pd.read_table("results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_TWSig.txt",sep=' ')
  chr_sig=sig_df.CHR.unique()

rule run_conditional:
  resources: mem_mb=70000 
  input:
    sig= "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_TWSig.txt",
    ss= "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_munged.sumstats.gz",
    fusion=rules.install_twas_plotter.output
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
      --save_loci \
      --ldsc F \
      --locus_win 500000"

if sig_file.is_file():
  rule conditional:
      input: expand("results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{chr}.joint_dropped.dat", chr=chr_sig)

######
# Process TWAS results
######

rule process_twas:
  input:
    "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt"
  output:
    "results/twas/twas_results/PGC_MDD3_twas_AllTissues_CLEAN.txt", "results/twas/twas_results/PGC_MDD3_twas_AllTissues_TWSig_CLEAN.txt", "results/twas/twas_results/PGC_MDD3_twas_panel_N.csv", "results/twas/twas_results/PGC_MDD3_twas.pos"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/process_twas.R"

######
# Create locus plots
######

rule run_plot_loci:
  input:
    twas_res= "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt",
    glist="resources/twas/glist-hg19",
    pos="results/twas/twas_results/PGC_MDD3_twas.pos", 
    fusion=rules.install_twas_plotter.output
  output:
    "results/twas/conditional/locus_plot_chr{chr}.log"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript resources/twas/TWAS-plotter/TWAS-locus-plotter.V1.0.r \
      --twas {input.twas_res} \
      --pos {input.pos} \
      --window 500000 \
      --gene_loc {input.glist} \
      --post_proc_prefix results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{wildcards.chr} | echo Done > results/twas/conditional/locus_plot_chr{wildcards.chr}.log"

if sig_file.is_file():
  rule plot_loci:
      input: expand("results/twas/conditional/locus_plot_chr{chr}.log", chr=chr_sig)

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

if sig_file.is_file():
  rule process_conditional:
    input: expand("results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{chr}.joint_dropped.dat", chr=chr_sig), "results/twas/twas_results/PGC_MDD3_twas_AllTissues_TWSig_CLEAN.txt"
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
  input: "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_munged.sumstats.gz"
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
    input: expand("results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv", chr=chr)

rule run_focus_gw:
  input: "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_munged.sumstats.gz"
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
    input: expand("results/twas/focus/GW_sig/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv", chr=chr)

# Process the FOCUS results
rule process_focus:
  input: expand("results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv", chr=chr)
  output:
    "results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.results.csv"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/process_focus.R"

########
# Create high confidence associations table
########

rule process_focus:
  input:
    fusion="results/twas/conditional/PGC_MDD3_TWAS_Conditional_table_novelty.csv", focus_res="results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.results.csv"
  output:
    "results/twas/PGC3_MDD_TWAS_HighConf_results.csv"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/create_high_conf.R"

########
# TWAS-GSEA
########

# Split twas results by tissue type
rule prep_for_gsea:
  input: "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt"
  output:
    "results/twas/twas_results/PGC_MDD3_twas_BLOOD_GW.txt","results/twas/twas_results/PGC_MDD3_twas_BRAIN_GW.txt","results/twas/twas_results/PGC_MDD3_twas_HPA_GW.txt","results/twas/twas_results/PGC_MDD3_twas_HPT_GW.txt","results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/prep_for_gsea.R"

# Run twas-gsea across results
rule run_twas_gsea:
  resources: mem_mb=30000, cpus=3
  input: ss="results/twas/twas_results/PGC_MDD3_twas_{set}_GW.txt", gmt="resources/twas/gsea/c2_and_c5.all.v7.4.entrez.gmt", twas_gsea=rules.install_twas_gsea.output
  output:
    "results/twas/gsea/PGC_MDD3_twas_{set}_GSEA.log"
  conda:
    "../envs/twas_gsea.yaml"
  shell:
    "Rscript {input.twas_gsea}/TWAS-GSEA.V1.2.R \
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
  input: ss="results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_munged.sumstats.gz", expr_scores=rules.download_mesc_score.output
  output: "results/twas/mesc/PGC_MDD3_TWAS.MESC.all.h2med"
  conda: "../envs/mesc.yaml"
  shell: "resources/twas/mesc/run_mesc.py --h2med {input.ss} --exp-chr {input.expr_scores}/All_Tissues --out results/twas/mesc/PGC_MDD3_TWAS.MESC"

# Estimate gene-set enrichment
rule est_h2med_sets:
  resources: mem_mb=20000 
  input: ss="results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp_munged.sumstats.gz", expr_score_sets=rules.download_mesc_score_sets.output
  output: "results/twas/mesc/PGC_MDD3_TWAS.MESC.sets.all.h2med"
  conda: "../envs/mesc.yaml"
  shell: "resources/twas/mesc/run_mesc.py --h2med {input.ss} --exp-chr {input.expr_score_sets}/All_Tissues --out results/twas/mesc/PGC_MDD3_TWAS.MESC.sets"
  
