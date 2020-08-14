import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# Copy summary statistics listed in config.yaml under sumstats
# with key FORMAT_COHORT.POP.hgNN.VERSION
rule stage_sumstats:
	input: lambda wildcards: config["sumstats"][wildcards.cohort]
	output: "resources/sumstats/{cohort}.gz"
	log: "logs/sumstats/stage/{cohort}.log"
	shell: "cp -v {input} {output} > {log}"

# Harmonize names of all summary statistics listed under sumstats in config.yaml
rule sumstats:
	input: expand("resources/sumstats/{sumstats}.gz", sumstats=config["sumstats"])
	
# For pre-formatted daner files
rule daner:
	input: "resources/sumstats/daner_{cohort}.gz"
	output: "results/sumstats/daner/daner_{cohort}.gz"
	log: "logs/sumstats/daner/daner_{cohort}.log"
	shell: "cp -v {input} {output} > {log}"
	
# Convert text sumstats to daner
rule text2daner:
	input: sumstats="resources/sumstats/text_mdd_{cohort}.{ancestries}.{build}.{version}.gz", sh="scripts/sumstats/{cohort}.sh"
	output: "results/sumstats/daner/daner_mdd_{cohort}.{ancestries}.{build}.{version}.gz"
	log: "logs/sumstats/daner/daner_mdd_{cohort}.{ancestries}.{build}.{version}.log"
	conda: "../envs/meta.yaml" 
	shell: "sh {input.sh} {input.sumstats} {output} {log}"
	
# for daner files on genome build hg19
rule hg19:
	input: "results/sumstats/daner/daner_mdd_{cohort}.{ancestries}.hg19.{version}.gz"
	output: "results/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{version}.gz"
	log: "logs/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{version}.log"
	shell: "cp -v {input} {output} > {log}"

# download hgIN to hgOUT chain
rule hg_chain:
    input: HTTP.remote("hgdownload.soe.ucsc.edu/goldenPath/hg{from}/liftOver/hg{from}ToHg{to}.over.chain.gz", keep_local=True)
	output: "resources/liftOver/hg{from}ToHg{to}.over.chain"
	run:
		 outputName = os.path.basename(input[0])
		 shell("gunzip -c {input} > {output}")
	
# liftover hg38 to hg19	
rule hg38to19:
	input: daner="results/sumstats/daner/daner_mdd_{cohort}.{ancestries}.hg38.{version}.gz", chain="resources/liftOver/hg38ToHg19.over.chain"
	output: "results/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{version}.gz"
	log: "logs/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{version}.log"
	conda: "../envs/meta.yaml" 
	script: "../scripts/liftover.R"

# align to imputation panel
rule align:
	input: daner="results/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.{build}.{version}.gz", ref="results/meta/reference_info"
	output: "results/sumstats/aligned/daner_mdd_{cohort}.{ancestries}.{build}.{version}.aligned.gz"
	log: "logs/sumstats/aligned/daner_mdd_{cohort}.{ancestries}.{build}.{version}.aligned.log"
	conda: "../envs/meta.yaml" 
	script: "../scripts/align.R"

# create reference into file linking to imputation panel
rule refdir:
	output: "results/meta/reference_info"
	log: "logs/meta/reference_info.log"
	shell: "cd results/meta; impute_dirsub --refdir {config[refdir]} --reference_info --outname meta"

# link sumstats files into meta-analysis directory
rule meta:
	input: "results/sumstats/aligned/{cohort}.gz"
	output: "results/meta/{cohort}.gz"
	log: "logs/meta/{cohort}.log"
	shell: "cp -v {input} {output} > {log}"

# Ricopili results dataset list for eur ancestries
rule dataset_eur:
	input: "results/meta/daner_mdd_MDD29.eur.hg19.0120a_rmUKBB.aligned.gz",
	 "results/meta/daner_mdd_23andMe.eur.hg19.v7_2.aligned.gz",
	 "results/meta/daner_mdd_deCODE.eur.hg19.160211.aligned.gz",
	 "results/meta/daner_mdd_GenScot.eur.hg19.1215a.aligned.gz",
	 "results/meta/daner_mdd_GERA.eur.hg19.0915a_mds5.aligned.gz",
	 "results/meta/daner_mdd_UKBB.eur.hg19.MD_glm.aligned.gz",
	 "results/meta/daner_mdd_iPSYCH.eur.hg19.170220.aligned.gz",
	 "results/meta/daner_mdd_FinnGen.eur.hg19.R5_18032020.aligned.gz",
	 "results/meta/daner_mdd_ALSPAC.eur.hg19.12082019.aligned.gz",
	 "results/meta/daner_mdd_Airwave.eur.hg19.0820.aligned.gz"
	output: "results/meta/dataset_full_eur_v{analysis}"
	log: "logs/meta/dataset_full_eur_v{analysis}.log"
	shell: "for daner in {input}; do echo $(basename $daner) >> {output}; done"
	
# Dataset list that exclude a particular cohort
rule dataset_nocCOHORT_eur:
	input: "results/meta/dataset_full_eur_v{analysis}"
	log: "logs/meta/dataset_no{cohort}_eur_v{analysis}"
	output: "results/meta/dataset_no{cohort}_eur_v{analysis}"
	shell: "cat {input} | grep --invert {wildcards.cohort} > {output}"
	

# Ricopili submission
rule postimp:
	input: dataset="results/meta/dataset_{cohorts}_{ancestries}_v{version}", ref="results/meta/reference_info"
	params:
		popname=lambda wildcards: wildcards.ancestries.upper(),
		dataset=lambda wildcards, input: os.path.basename(input.dataset)
	output: touch("results/meta/{cohorts}_{ancestries}_v{version}.done")
	log: "logs/meta/pgc_mdd_meta_{cohorts}_{ancestries}_hg19_v{version}.postimp_navi.log"
	shell: "cd results/meta; postimp_navi --result {params.dataset} --popname {params.popname} --nolahunt --out pgc_mdd_{cohorts}_{wildcards.ancestries}_hg19_v{wildcards.version}"

# current European ancestries analysis
# analysis version format: v3.[PGC Cohorts Count].[Other Cohorts Count]_YYYY-MM-DD
rule postimp_eur:
	input: "results/meta/full_eur_v3.29.09.done"
	

# distribute results

# glob of all files in the distribution directory
distribution_full_gz, = glob_wildcards("results/meta/distribution/pgc_mdd_full_eur_hg19_v3.29.08/{file}.gz")
distribution_full_xls, = glob_wildcards("results/meta/distribution/pgc_mdd_full_eur_hg19_v3.29.08/{file}.xls")
distribution_full_pdf, = glob_wildcards("results/meta/distribution/pgc_mdd_full_eur_hg19_v3.29.08/{file}.pdf")
# check if glob doesn't return anything, and assign nonsense values
# this allows the expand() statement in the dbox_full fule to not fail even though
# the rule won't actually be run
distribution_full_gz = distribution_full_gz if distribution_full_gz else ['spurious']
distribution_full_xls = distribution_full_xls if distribution_full_xls else ['spurious']
distribution_full_pdf = distribution_full_pdf if distribution_full_pdf else ['spurious']

rule distribute_full:
	input: "results/meta/distribution/pgc_mdd_full_eur_hg19_v3.29.08/{file}"
	output: DBox.remote("distribution/pgc_mdd_full_eur_hg19_v3.29.08/{file}")
	shell: "cp {input} {output}"

# list all files to be uploaded to Dropbox
rule dbox_full:
	input: DBox.remote(expand("distribution/pgc_mdd_full_eur_hg19_v3.29.08/{file}.gz", file=distribution_full_gz)), \
	       DBox.remote(expand("distribution/pgc_mdd_full_eur_hg19_v3.29.08/{file}.xls", file=distribution_full_xls)), \
		   DBox.remote(expand("distribution/pgc_mdd_full_eur_hg19_v3.29.08/{file}.pdf", file=distribution_full_pdf))

# Download full sumstats for downstream analysis
rule redistribute_full:
	input: DBox.remote("distribution/{analysis}/daner_{analysis}.gz")
	output: "results/distribution/daner_{analysis}.gz"
	shell: "cp {input} {output}"

rule downstream_full:
	input: "results/distribution/daner_pgc_mdd_full_eur_hg19_v3.29.08.gz"
