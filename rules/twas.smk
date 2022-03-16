##########
# Run twas analysis
##########

####
# Install required software not within conda
####

# Install fusion
rule install_fusion:
  output:
    directory("resources/twas/fusion/")
  conda:
    "../envs/twas.yaml"
  shell:
    "git clone git@github.com:gusevlab/fusion_twas.git {output}"

# Download plink2R
rule download_plink2R:
  input:
    rules.install_fusion.output
  output:
    "resources/twas/plink2R/plink2R-master/data.bed"
  conda:
    "../envs/twas.yaml"
  shell:
    "wget -O resources/twas/plink2R/master.zip https://github.com/gabraham/plink2R/archive/master.zip; unzip resources/twas/plink2R/master.zip -d resources/twas/plink2R"

# Install plink2R
rule install_plink2R:
  input: 
    rules.download_plink2R.output
  output:
    touch("resources/twas/install_plink2R")
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript -e 'install.packages(\"resources/twas/plink2R/plink2R-master/plink2R/\",repos=NULL)'"

# Install focus
rule install_focus:
  conda:
    "../envs/twas.yaml"
  output: 
    touch("resources/twas/pyfocus")
  shell:
    "pip install pyfocus==0.6.10 --user"

# Install TWAS-plotter
rule install_twas_plotter:
  output:
    directory("resources/twas/TWAS-plotter/")
  conda:
    "../envs/twas.yaml"
  shell:
    "git clone git@github.com:opain/TWAS-plotter.git {output}"

# Install SNP-weights pipeline repo
rule install_snp_weight_pipe:
  output:
    directory("resources/twas/Calculating-FUSION-TWAS-weights-pipeline/")
  conda:
    "../envs/twas.yaml"
  shell:
    "git clone git@github.com:opain/Calculating-FUSION-TWAS-weights-pipeline.git {output}"

# Install MESC
rule install_mesc:
  output:
    directory("resources/twas/mesc/")
  conda:
    "../envs/twas.yaml"
  shell:
    "git clone git@github.com:douglasyao/mesc.git {output}"

# Install ggvenn
rule install_ggvenn:
  output:
    touch("resources/software/install_ggvenn.done")
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript -e 'devtools::install_github(\"yanlinlin82/ggvenn\")'"

####
# Download required data
####

# Create a list of FUSION weights
weights=["Adrenal_Gland","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Substantia_nigra","Pituitary","Thyroid","Whole_Blood","CMC.BRAIN.RNASEQ","CMC.BRAIN.RNASEQ_SPLICING","NTR.BLOOD.RNAARR","YFS.BLOOD.RNAARR"]

# Create list of FUSION SNP-weight sets from GTEx
gtex_weights=["Adrenal_Gland","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Substantia_nigra","Pituitary","Thyroid","Whole_Blood"]

# Create list of FUSION SNP-weight sets not from GTEx
non_gtex_weights=["CMC.BRAIN.RNASEQ","CMC.BRAIN.RNASEQ_SPLICING","NTR.BLOOD.RNAARR","YFS.BLOOD.RNAARR"]

# Download FUSION GTEx SNP-weights
rule download_gtex_weights:
  output:
    "resources/twas/fusion_data/{weight}/{weight}.pos"
  conda:
    "../envs/twas.yaml"
  shell:
    "wget -O resources/twas/fusion_data/{wildcards.weight}.P01.tar.bz2 http://gusevlab.org/projects/fusion/weights/GTEx.{wildcards.weight}.P01.tar.bz2; tar xvjf resources/twas/fusion_data/{wildcards.weight}.P01.tar.bz2 -C resources/twas/fusion_data/{wildcards.weight}; rm resources/twas/fusion_data/{wildcards.weight}.P01.tar.bz2; mv resources/twas/fusion_data/{wildcards.weight}/{wildcards.weight}.P01.pos resources/twas/fusion_data/{wildcards.weight}/{wildcards.weight}.pos"

rule gtex_weights:
    input: expand("resources/twas/fusion_data/{weight}/{weight}.pos", weight=gtex_weights)

# Download FUSION non-GTEx SNP-weights
rule download_non_gtex_weights:
  output:
    "resources/twas/fusion_data/{weight}/{weight}.profile"
  conda:
    "../envs/twas.yaml"
  shell:
    "wget -O resources/twas/fusion_data/{wildcards.weight}.tar.bz2 https://data.broadinstitute.org/alkesgroup/FUSION/WGT/{wildcards.weight}.tar.bz2; tar xvjf resources/twas/fusion_data/{wildcards.weight}.tar.bz2 -C resources/twas/fusion_data/{wildcards.weight}; rm resources/twas/fusion_data/{wildcards.weight}.tar.bz2"

rule non_gtex_weights:
  input: expand("resources/twas/fusion_data/{weight}/{weight}.profile", weight=non_gtex_weights)

# Insert N into non-GTEX SNP-weights
rule insert_n_nongtex:
  input:
    rules.non_gtex_weights.input
  output:
    touch("resources/twas/insert_n_nongtex.done")
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/insert_n_nongtex.R"

# Download PsychENCODE SNP-weights
rule download_psychENCODE_weights:
  output:
    directory("resources/twas/psychencode_data/SNP-weights/PEC_TWAS_weights")
  conda:
    "../envs/twas.yaml"
  shell:
    "wget -O resources/twas/psychencode_data/PEC_TWAS_weights.tar.gz http://resource.psychencode.org/Datasets/Derived/PEC_TWAS_weights.tar.gz; mkdir -p resources/twas/psychencode_data/SNP-weights/PEC_TWAS_weights; tar xvzf resources/twas/psychencode_data/PEC_TWAS_weights.tar.gz -C resources/twas/psychencode_data/SNP-weights/PEC_TWAS_weights; rm resources/twas/psychencode_data/PEC_TWAS_weights.tar.gz"

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
# Download and format LD reference data
###

# Create list of chromosome numbers
chr=range(1, 23)

rule prep_1kg:
  resources:
    mem_mb=20000
  output:
    expand("resources/twas/1kg/EUR/EUR_phase3.MAF_001.chr{chr}.bed",chr=chr)
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/1kg_prep.R"

###
# Format PsychENCODE SNP-weights
###

rule format_psychencode:
  input:
    psychencode_data=rules.download_psychENCODE_weights.output, 
    weights_pipe=rules.install_snp_weight_pipe.output
  output:
    touch("resources/twas/psychencode_data/format_psychencode.done")
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/format_psychENCODE.R"

###
# Create database for FOCUS
###

rule focus_db:
  input:
    gtex_weights=rules.gtex_weights.input,
    non_gtex_weights=rules.insert_n_nongtex.output,
    psychencode_data=rules.format_psychencode.output,
    focus=rules.install_focus.output
  output:
    "resources/twas/focus_db/MDD_TWAS.db"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/focus_db_prep.R"

###
# Munge sumstats
###

# Format sumstats file so FOCUS can read it
rule pre_munge:
  input:
    "results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.10.neff.gz"
  output:
    "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.10.neff_premunged.gz"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/pre_munge.R"
    
# munge sumstats using FOCUS munge function
rule focus_munge:
  input:
    premunged=rules.pre_munge.output,
    focus=rules.install_focus.output
  output:
    "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.10.neff_munged.sumstats.gz"
  conda:
    "../envs/twas.yaml"
  shell:
    "focus munge {input.premunged} --output results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.10.neff_munged"

###
# Run TWAS
###

# Calculate median effective sample size
# FUSION requires this parameter to be specified despite having the N column in the sumstats
rule retrieve_Neff:
  input:
    rules.focus_munge.output
  output:
    "results/twas/median_Neff.txt"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/median_neff.R --munged {input} --out {output}"

# run twas
rule run_twas:
  resources:
    mem_mb=20000
  input:
    sumstats=rules.focus_munge.output,
    neff_txt="results/twas/median_Neff.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    gtex_weights=rules.gtex_weights.input,
    non_gtex_weights=rules.insert_n_nongtex.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "results/twas/PGC_MDD3_twas_{weight}_chr{chr}"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Neff=$(cat results/twas/median_Neff.txt); Rscript resources/twas/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights resources/twas/fusion_data/{wildcards.weight}/{wildcards.weight}.pos \
    --weights_dir resources/twas/fusion_data/{wildcards.weight} \
    --ref_ld_chr resources/twas/1kg/EUR/EUR_phase3.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 1e-3 \
    --GWASN ${{Neff}}"

rule fusion_twas:
    input: expand("results/twas/PGC_MDD3_twas_{weight}_chr{chr}", weight=weights, chr=chr)

# Run twas using psychENCODE SNP-weights
rule run_psychencode_twas:
  resources: mem_mb=20000 
  input:
    sumstats=rules.focus_munge.output,
    neff_txt="results/twas/median_Neff.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    format_psychencode=rules.format_psychencode.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "results/twas/psychencode/PGC_MDD3_twas_psychencode_chr{chr}"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Neff=$(cat results/twas/median_Neff.txt); Rscript resources/twas/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights resources/twas/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights.pos \
    --weights_dir resources/twas/psychencode_data/PEC_TWAS_weights \
    --ref_ld_chr resources/twas/1kg/EUR/EUR_phase3.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 1e-3 \
    --GWASN ${{Neff}}"

rule psychencode_twas:
    input: expand("results/twas/psychencode/PGC_MDD3_twas_psychencode_chr{chr}", chr=chr)

###
# Combine TWAS results
###

rule combine_twas_res:
  input:
    rules.fusion_twas.input,
    rules.psychencode_twas.input
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
    results=rules.combine_twas_res.output,
    twas_plotter=rules.install_twas_plotter.output
  output:
    "results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW_Manhattan.png"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript resources/twas/TWAS-plotter/TWAS-plotter.V1.0.r \
    --twas results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt \
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
rule run_conditional:
  resources: mem_mb=70000 
  input:
    sig=rules.combine_twas_res.output,
    ss=rules.focus_munge.output,
    glist=rules.download_glist.output,
    fusion=rules.install_fusion.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{chr}.joint_dropped.dat"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript /scratch/groups/biomarkers-brc-mh/TWAS_resource/FUSION/fusion_twas-master/FUSION.post_process.R \
      --input {input.sig} \
      --sumstats {input.ss} \
      --report \
      --ref_ld_chr resources/twas/1kg/EUR/EUR_phase3.MAF_001.chr \
      --out results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{wildcards.chr} \
      --chr {wildcards.chr} \
      --save_loci \
      --ldsc F \
      --locus_win 500000"

rule conditional:
      input: expand("results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{chr}.joint_dropped.dat", chr=chr)

######
# Process TWAS results
######

rule process_twas:
  input:
    rules.combine_twas_res.output
  output:
    "results/twas/twas_results/PGC_MDD3_twas_AllTissues_CLEAN.txt", 
    "results/twas/twas_results/PGC_MDD3_twas_AllTissues_TWSig_CLEAN.txt", 
    "results/twas/twas_results/PGC_MDD3_twas_panel_N.csv", 
    "results/twas/twas_results/PGC_MDD3_twas.pos"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/process_twas.R"

######
# Create locus plots
######

rule run_plot_loci:
  input:
    twas_res=rules.conditional.input,
    glist=rules.download_glist.output,
    process_twas=rules.process_twas.output,
    fusion=rules.install_twas_plotter.output
  output:
    touch("results/twas/conditional/locus_plot_chr{chr}.done")
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript resources/twas/TWAS-plotter/TWAS-locus-plotter.V1.0.r \
      --twas results/twas/twas_results/PGC_MDD3_twas_AllTissues_GW.txt \
      --pos results/twas/twas_results/PGC_MDD3_twas.pos \
      --window 500000 \
      --gene_loc {input.glist} \
      --post_proc_prefix results/twas/conditional/PGC_MDD3_TWAS_conditional_chr{wildcards.chr}"

rule plot_loci:
    input: expand("results/twas/conditional/locus_plot_chr{chr}.done", chr=chr)

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
  input: 
    rules.conditional.input,
    rules.process_twas.output
  output:
    "results/twas/conditional/PGC_MDD3_TWAS_Conditional_table_novelty.csv"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/process_conditional.R"

########
# Finemap TWAS associations
########

# Note I was having a conda compatability error
# So I used a local installation of focus in the end
rule run_focus:
  input:
    sumstats=rules.focus_munge.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv"
#  conda: 
#    "../envs/twas.yaml"
  shell:
    "/users/k1806347/brc_scratch/Software/focus-master/bin/focus finemap {input.sumstats} \
    	resources/twas/1kg/EUR/EUR_phase3.MAF_001.chr{wildcards.chr} \
    	/scratch/groups/biomarkers-brc-mh/TWAS_resource/FOCUS/MDD_TWAS_db/MDD_TWAS.db \
    	--chr {wildcards.chr} \
    	--p-threshold 1e-4 \
    	--out results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{wildcards.chr}"

# Check FOCUS completed
rule check_focus:
  input:
    "results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv"
  output:
    "results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.complete"
  conda: 
    "../envs/twas.yaml"
  shell:
    "if [ $(grep -c ERROR results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{wildcards.chr}.log) -eq 0 ]; then touch results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{wildcards.chr}.focus.complete; else rm results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{wildcards.chr}.focus.tsv; fi"

rule focus:
    input: expand("results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.complete", chr=chr)

# Run FOCUS using only genome-wide significance threshold out of interest
rule run_focus_gw:
  input: 
    sumstats=rules.focus_munge.output,
    prep_1kg=rules.prep_1kg.output,
    focus_db=rules.focus_db.output
  output: 
    "results/twas/focus/GW_sig/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv"
  conda: 
    "../envs/twas.yaml"
  shell:
    "focus finemap {input.sumstats} \
    	resources/twas/1kg/EUR/EUR_phase3.MAF_001.chr{wildcards.chr} \
    	resources/twas/focus_db/MDD_TWAS.db \
    	--chr {wildcards.chr} \
    	--out results/twas/focus/GW_sig/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{wildcards.chr}"

rule focus_gw:
    input: expand("results/twas/focus/GW_sig/PGC_MDD3_TWAS.FOCUS.MDD_TWAS_db.chr{chr}.focus.tsv", chr=chr)

# Process the FOCUS results
rule process_focus:
  input: 
    rules.focus.input,
    rules.process_twas.output
  output:
    "results/twas/focus/p_1e-4/PGC_MDD3_TWAS.FOCUS.results.csv"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/process_focus.R"

########
# Create high confidence associations table
########

rule process_high_conf:
  input:
    fusion=rules.process_conditional.output,
    focus_res=rules.process_focus.output
  output:
    "results/twas/PGC3_MDD_TWAS_HighConf_results.csv"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/create_high_conf.R"

########
# Run MESC
########

# Estimate h2med
rule est_h2med:
  input: 
    ss=rules.focus_munge.output,
    expr_scores=rules.download_mesc_score.output
  output: 
    "results/twas/mesc/PGC_MDD3_TWAS.MESC.all.h2med"
  conda: 
    "../envs/mesc.yaml"
  shell: 
    "resources/twas/mesc/run_mesc.py --h2med {input.ss} --exp-chr {input.expr_scores}/All_Tissues --out results/twas/mesc/PGC_MDD3_TWAS.MESC"

# Estimate gene-set enrichment
rule est_h2med_sets:
  resources: 
    mem_mb=20000 
  input: 
    ss=rules.focus_munge.output,
    expr_score_sets=rules.download_mesc_score_sets.output
  output: 
    "results/twas/mesc/PGC_MDD3_TWAS.MESC.sets.all.h2med"
  conda: 
    "../envs/mesc.yaml"
  shell: 
    "resources/twas/mesc/run_mesc.py --h2med {input.ss} --exp-chr {input.expr_score_sets}/All_Tissues --out results/twas/mesc/PGC_MDD3_TWAS.MESC.sets"

########
# SMR using eQTLGen and MetaBrain
########

# Download SMR
rule download_smr:
  output: 
    directory("resources/twas/smr")
  conda: 
    "../envs/mesc.yaml"
  shell: 
    "mkdir resources/twas/smr; wget -c --tries=0 --read-timeout=20 --no-check-certificate -O resources/twas/smr/smr_Linux.zip https://cnsgenomics.com/software/smr/download/smr_Linux.zip; unzip resources/twas/smr/smr_Linux.zip -d resources/twas/smr/; rm resources/twas/smr/smr_Linux.zip"
  
# Download eQTLGen data in SMR format
rule download_eQTLGen:
  output: 
    directory("resources/twas/eQTLGen")
  conda: 
    "../envs/mesc.yaml"
  shell: 
    "mkdir resources/twas/eQTLGen; wget -O resources/twas/eQTLGen/cis-eQTL-SMR_20191212.tar.gz https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz; tar -xvzf resources/twas/eQTLGen/cis-eQTL-SMR_20191212.tar.gz -C resources/twas/eQTLGen/; rm resources/twas/eQTLGen/cis-eQTL-SMR_20191212.tar.gz; gunzip resources/twas/eQTLGen/*"

##
# Download MetaBrain data in SMR format
##

# Basalganglia
rule download_MetaBrain_Basalganglia:
  output: 
    directory("resources/twas/MetaBrain/Basalganglia")
  conda: 
    "../envs/mesc.yaml"
  shell: 
    "mkdir -p resources/twas/MetaBrain/Basalganglia; wget -O resources/twas/MetaBrain/Basalganglia/2020-05-26-Basalganglia-EUR-smr.zip https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-Basalganglia-EUR/2020-05-26-Basalganglia-EUR-smr.zip; unzip -d resources/twas/MetaBrain/Basalganglia/ resources/twas/MetaBrain/Basalganglia/2020-05-26-Basalganglia-EUR-smr.zip; rm resources/twas/MetaBrain/Basalganglia/2020-05-26-Basalganglia-EUR-smr.zip"

# Cerebellum
rule download_MetaBrain_Cerebellum:
  output: 
    directory("resources/twas/MetaBrain/Cerebellum")
  conda: 
    "../envs/mesc.yaml"
  shell: 
    "mkdir -p resources/twas/MetaBrain/Cerebellum; wget -O resources/twas/MetaBrain/Cerebellum/2020-05-26-Cerebellum-EUR-smr.zip https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-Cerebellum-EUR/2020-05-26-Cerebellum-EUR-smr.zip; unzip -d resources/twas/MetaBrain/Cerebellum/ resources/twas/MetaBrain/Cerebellum/2020-05-26-Cerebellum-EUR-smr.zip; rm resources/twas/MetaBrain/Cerebellum/2020-05-26-Cerebellum-EUR-smr.zip"

# Cortex
rule download_MetaBrain_Cortex:
  output: 
    directory("resources/twas/MetaBrain/Cortex")
  conda: 
    "../envs/mesc.yaml"
  shell: 
    "mkdir -p resources/twas/MetaBrain/Cortex; wget -O resources/twas/MetaBrain/Cortex/2020-05-26-Cortex-EUR-smr.zip https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-Cortex-EUR/2020-05-26-Cortex-EUR-smr.zip; unzip -d resources/twas/MetaBrain/Cortex/ resources/twas/MetaBrain/Cortex/2020-05-26-Cortex-EUR-smr.zip; rm resources/twas/MetaBrain/Cortex/2020-05-26-Cortex-EUR-smr.zip"

# Hippocampus
rule download_MetaBrain_Hippocampus:
  output: 
    directory("resources/twas/MetaBrain/Hippocampus")
  conda: 
    "../envs/mesc.yaml"
  shell: 
    "mkdir -p resources/twas/MetaBrain/Hippocampus; wget -O resources/twas/MetaBrain/Hippocampus/2020-05-26-Hippocampus-EUR-smr.zip https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-Hippocampus-EUR/2020-05-26-Hippocampus-EUR-smr.zip; unzip -d resources/twas/MetaBrain/Hippocampus/ resources/twas/MetaBrain/Hippocampus/2020-05-26-Hippocampus-EUR-smr.zip; rm resources/twas/MetaBrain/Hippocampus/2020-05-26-Hippocampus-EUR-smr.zip"

# Spinalcord
rule download_MetaBrain_Spinalcord:
  output: 
    directory("resources/twas/MetaBrain/Spinalcord")
  conda: 
    "../envs/mesc.yaml"
  shell: 
    "mkdir -p resources/twas/MetaBrain/Spinalcord; wget -O resources/twas/MetaBrain/Spinalcord/2020-05-26-Spinalcord-EUR-smr.zip https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-Spinalcord-EUR/2020-05-26-Spinalcord-EUR-smr.zip; unzip -d resources/twas/MetaBrain/Spinalcord/ resources/twas/MetaBrain/Spinalcord/2020-05-26-Spinalcord-EUR-smr.zip; rm resources/twas/MetaBrain/Spinalcord/2020-05-26-Spinalcord-EUR-smr.zip"

rule download_MetaBrain_all:
  input: 
    rules.download_MetaBrain_Basalganglia.output,
    rules.download_MetaBrain_Cerebellum.output,
    rules.download_MetaBrain_Cortex.output,
    rules.download_MetaBrain_Hippocampus.output,
    rules.download_MetaBrain_Spinalcord.output
  output:
    touch('resources/twas/MetaBrain_download.out')

# Update variant IDs in MetaBrain SMR files
rule format_metabrain_esi:
  output: 
    touch("resources/twas/MetaBrain/format_MetaBrain_esi.out")
  conda: 
    "../envs/twas.yaml"
  shell: 
    "Rscript scripts/twas/format_metabrain_esi.R"

# Convert sumstats to COJO format
rule daner_to_cojo:
  output:
    "results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.10.neff_COJO.txt"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/daner_to_cojo.R"

# Run SMR with eQTLGen
rule smr_analysis_eQTLGen:
  resources:
    mem_mb=15000 
  input:
    rules.prep_1kg.output,
    rules.daner_to_cojo.output,
    rules.download_eQTLGen.output,
    rules.download_smr.output
  output:
    "results/twas/eqtlgen_smr/eqtlgen_smr_res_chr{chr}.smr"
  conda: 
    "../envs/twas.yaml"
  shell:
    "/users/k1806347/brc_scratch/Software/smr_Linux \
      --bfile resources/twas/1kg/all_phase3.chr{wildcards.chr} \
      --keep resources/twas/1kg/super_pop_keep_files/EUR_samples.keep \
      --gwas-summary results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.10.neff_COJO.txt \
      --beqtl-summary resources/twas/eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
      --out results/twas/eqtlgen_smr/eqtlgen_smr_res_chr{wildcards.chr} \
      --thread-num 1"

rule run_smr_analysis_eQTLGen:
  input: expand('results/twas/eqtlgen_smr/eqtlgen_smr_res_chr{chr}.smr', chr=chr)
  output: touch('results/twas/eqtlgen_smr/smr_eQTLGen_complete.out')

# Run SMR with MetaBrain
rule smr_analysis_MetaBrain:
  resources:
    mem_mb=15000 
  input:
    rules.prep_1kg.output,
    rules.daner_to_cojo.output,
    rules.format_metabrain_esi.output,
    rules.download_smr.output
  output:
    "results/twas/metabrain_smr/metabrain_{tissue}_smr_res_chr{chr}.smr"
  conda: 
    "../envs/twas.yaml"
  shell:
    "/users/k1806347/brc_scratch/Software/smr_Linux \
      --bfile resources/twas/1kg/all_phase3.chr{wildcards.chr} \
      --keep resources/twas/1kg/super_pop_keep_files/EUR_samples.keep \
      --gwas-summary results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.10.neff_COJO.txt \
      --beqtl-summary resources/twas/MetaBrain/{wildcards.tissue}/2020-05-26-{wildcards.tissue}-EUR-{wildcards.chr}-SMR-besd \
      --out results/twas/metabrain_smr/metabrain_{wildcards.tissue}_smr_res_chr{wildcards.chr} \
      --thread-num 1"

rule run_smr_analysis_MetaBrain:
  input: expand('results/twas/metabrain_smr/metabrain_{tissue}_smr_res_chr{chr}.smr', chr=chr, tissue=["Basalganglia","Cerebellum","Cortex","Hippocampus","Spinalcord"])
  output: touch('results/twas/metabrain_smr/smr_metabrain_complete.out')

# Process SMR results
rule process_smr:
  input:
    rules.run_smr_analysis_eQTLGen.output,
    rules.run_smr_analysis_MetaBrain.output
  output:
    "results/twas/eqtlgen_smr/eqtlgen_smr_res_GW_withIDs.csv"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/twas/process_smr.R"

########
# Create report of the results
########

rule create_report:
  input:
    rules.make_manhattan.output,
    rules.process_coloc.output,
    rules.process_high_conf.output,
    rules.plot_loci.input,
    rules.est_h2med.output,
    rules.est_h2med_sets.output,
    rules.process_smr.output,
    rules.install_ggvenn.output
  output:
    "docs/twas_report.html"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript -e \"rmarkdown::render(\'scripts/twas/twas_report.Rmd\', output_file = \'../../docs/twas_report.html\')\""
