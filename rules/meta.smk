# Conduct meta-analysis in Ricopili


#################
#               #
# Meta analyses #
#               #
#################

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
	shell: "cat {input} | grep --invert daner_mdd_{wildcards.cohort} > {output}"
	

# Ricopili submission
rule postimp:
	input: dataset="results/meta/dataset_{cohorts}_{ancestries}_v{version}", ref="results/meta/reference_info"
	params:
		popname=lambda wildcards: wildcards.ancestries.upper(),
		dataset=lambda wildcards, input: os.path.basename(input.dataset)
	output: touch("results/meta/{cohorts}_{ancestries}_v{version}.done")
	log: "logs/meta/pgc_mdd_meta_{cohorts}_{ancestries}_hg19_v{version}.postimp_navi.log"
	shell: "cd results/meta; postimp_navi --result {params.dataset} --popname {params.popname} --nolahunt --out pgc_mdd_{wildcards.cohorts}_{wildcards.ancestries}_hg19_v{wildcards.version}"

# current European ancestries analysis
# analysis version format: v3.[PGC Cohorts Count].[Other Cohorts Count]
analysis_version = ["3.29.09"]
rule postimp_eur:
	input: expand("results/meta/full_eur_v{version}.done", version=analysis_version)
	
# cohort sets for analysts
cohorts_analyst = ["full", "noUKBB", "no23andMe", "noALSPAC"]

# cohort sets for public
cohorts_public = ["no23andMe"]



###################################
#                                 #
# Summary statistics distribution #
#                                 #
###################################


# Distribute results
# extensions and prefixes of Ricopili distribution output files
# daner_pgc_mdd_full_eur_hg19_v{version}.EXT
distribution_daner_ext = ["gz", "gz.ldsc.sumstats.gz", "gz.p3.gz", "gz.p4.clump.areator.sorted.1mhc", "gz.p4.clump.areator.sorted.1mhc.pvsorted", "gz.p4.clump.areator.sorted.1mhc.pvsorted.regs.txt", "gz.p4.clump.areator.sorted.1mhc.summary", "gz.p4.clump.areator.sorted.1mhc.xls", "het.gz.p4.clump.areator.sorted.1mhc", "het.gz.p4.clump.areator.sorted.1mhc.xls"]

# PREFIX.pgc_mdd_full_eur_hg19_v{version}.pdf
distribution_pdf_prefix = ["areas.fo", "areas", "manhattan.nog2", "manhattan.nog", "manhattan.v2", "qq"]

# PREFIX.pgc_mdd_full_eur_hg19_v{version}.het.pdf
distribution_het_pdf_prefix = ["manhattan.v2", "qq"]

# basic.pgc_mdd_full_eur_hg19_v{version}.EXT
distribution_basic_ext = ["num.xls"]

# Distribute meta analysis files
rule distribute_meta:
	input: "results/meta/distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{file}"
	output: DBox_dist.remote("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{file}")
	shell: "cp {input} {output}"

##
## Distribution for PGC analysts
##

# list all files to be uploaded to Dropbox
rule DBox_dist_analyst:
	input: DBox_dist.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/daner_pgc_mdd_{cohorts}_eur_hg19_v{version}.{ext}", version=analysis_version, cohorts=cohorts_analyst, ext=distribution_daner_ext)), DBox_dist.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{prefix}.pgc_mdd_{cohorts}_eur_hg19_v{version}.pdf", version=analysis_version, cohorts=cohorts_analyst, prefix=distribution_pdf_prefix)), DBox_dist.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{prefix}.pgc_mdd_{cohorts}_eur_hg19_v{version}.het.pdf", version=analysis_version, cohorts=cohorts_analyst, prefix=distribution_het_pdf_prefix)), DBox_dist.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/basic.pgc_mdd_{cohorts}_eur_hg19_v{version}.num.xls", version=analysis_version, cohorts=cohorts_analyst, ext=distribution_basic_ext))

# Download daner sumstats for downstream analysis
rule redistribute_daner:
	input: DBox_dist.remote("distribution/{analysis}/daner_{analysis}.gz")
	output: "results/distribution/daner_{analysis}.gz"
	shell: "cp {input} {output}"

rule downstream_full:
	input: expand("results/distribution/daner_pgc_mdd_full_eur_hg19_v{version}.gz", version=analysis_version)

rule downstream_noUKBB:
	input: expand("results/distribution/daner_pgc_mdd_noUKBB_eur_hg19_v{version}.gz", version=analysis_version)

rule downstream_no23andMe:
	input: expand("results/distribution/daner_pgc_mdd_no23andMe_eur_hg19_v{version}.gz", version=analysis_version)

# Download tables and figures
rule redistribute_figtabs:
	input: DBox_dist.remote("distribution/{analysis}_v{version}/{prefix}.{analysis}_v{version}.{ext}")
	output: "results/distribution/{prefix}.{analysis}_v{version}.{ext}"
	shell: "cp {input} {output}"

# download most recent manhattan plot
rule manhattan_full:
	input: expand("results/distribution/manhattan.nog2.pgc_mdd_full_eur_hg19_v{version}.pdf", version=analysis_version)


##
## Distribution for public
##

ruleorder: distribute_public > distribute_meta

rule distribute_public:
	input: "results/meta/distribution/pgc_mdd_no23andMe_eur_hg19_v{version}/{file}"
	output: DBox_dist_public.remote("distribution/pgc_mdd_no23andMe_eur_hg19_v{version}/{file}")
	shell: "cp {input} {output}"

# list all files to be uploaded to Dropbox
rule DBox_dist_public:
	input: DBox_dist_public.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/daner_pgc_mdd_{cohorts}_eur_hg19_v{version}.{ext}", version=analysis_version, cohorts=cohorts_public, ext=distribution_daner_ext)), DBox_dist_public.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{prefix}.pgc_mdd_{cohorts}_eur_hg19_v{version}.pdf", version=analysis_version, cohorts=cohorts_public, prefix=distribution_pdf_prefix)), DBox_dist_public.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{prefix}.pgc_mdd_{cohorts}_eur_hg19_v{version}.het.pdf", version=analysis_version, cohorts=cohorts_public, prefix=distribution_het_pdf_prefix)), DBox_dist_public.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/basic.pgc_mdd_{cohorts}_eur_hg19_v{version}.num.xls", version=analysis_version, cohorts=cohorts_public, ext=distribution_basic_ext))


