import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# Copy summary statistics listed in config.yaml under sumstats
# with key FORMAT_COHORT.POP.hgNN.VERSION
rule stage_sumstats:
	input: lambda wildcards: config["sumstats"][wildcards.cohort]
	output: "resources/sumstats/{cohort}.gz"
	shell: "cp {input} {output}"

# Harmonize names of all summary statistics listed under sumstats in config.yaml
rule sumstats:
	input: expand("resources/sumstats/{sumstats}.gz", sumstats=config["sumstats"])
	
# For pre-formatted daner files
rule daner:
	input: "resources/sumstats/daner_{cohort}.gz"
	output: "results/sumstats/daner/daner_{cohort}.gz"
	shell: "cp {input} {output}"
	
# Convert text sumstats to daner
rule text2daner:
	input: sumstats="resources/sumstats/text_mdd_{cohort}.{ancestries}.{build}.{version}.gz", sh="scripts/sumstats/{cohort}.sh"
	output: "results/sumstats/daner/daner_mdd_{cohort}.{ancestries}.{build}.{version}.gz"
	shell: "sh {input.sh} {input.sumstats} {output}"
	
# for daner files on genome build hg19
rule hg19:
	input: "results/sumstats/daner/daner_mdd_{cohort}.{ancestries}.hg19.{version}.gz"
	output: "results/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{version}.gz"
	shell: "cp {input} {output}"

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
	script: "../scripts/liftover.R"

# align to imputation panel
rule align:
	input: daner="results/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.{build}.{version}.gz", ref="results/meta/reference_info"
	output: "results/sumstats/aligned/daner_mdd_{cohort}.{ancestries}.{build}.{version}.aligned.gz"
	script: "../scripts/align.R"

# create reference into file linking to imputation panel
rule refdir:
	output: "results/meta/reference_info"
	shell: "cd results/meta; impute_dirsub --refdir {config[refdir]} --reference_info --outname meta"

# link sumstats files into meta-analysis directory
rule meta:
	input: "results/sumstats/aligned/{cohort}.gz"
	output: "results/meta/{cohort}.gz"
	shell: "cp {input} {output}"

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
	 "results/meta/daner_mdd_ALSPAC.eur.hg19.12082019.aligned.gz"
	output: "results/meta/dataset_eur_v{analysis}"
	shell: "for daner in {input}; do echo $(basename $daner) >> {output}; done"

# Ricopili submission
rule postimp:
	input: dataset="results/meta/dataset_{ancestries}_v{analysis}", ref="results/meta/reference_info"
	params:
		popname=lambda wildcards: wildcards.ancestries.upper(),
		dataset=lambda wildcards, input: os.path.basename(input.dataset)
	output: touch("results/meta/{ancestries}_v{analysis}.done")
	shell: "cd results/meta; postimp_navi --result {params.dataset} --popname {params.popname} --nolahunt --out pgc_mdd_meta_{wildcards.ancestries}_hg19_v{wildcards.analysis}"

# current European ancestries analysis
# analysis version format: v3.[PGC Cohorts Count].[Other Cohorts Count]_YYYY-MM-DD
rule postimp_eur:
	input: "results/meta/eur_v3.29.08_2020-07-29.done"
	