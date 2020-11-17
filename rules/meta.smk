# Conduct meta-analysis in Ricopili


#################
#               #
# Meta analyses #
#               #
#################

# Copy summary statistics listed in config.yaml under sumstats
# with key FORMAT_COHORT.POP.hgNN.RELEASE
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
	input: sumstats="resources/sumstats/text_mdd_{cohort}.{ancestries}.{build}.{release}.gz", sh="scripts/sumstats/{cohort}.sh"
	output: "results/sumstats/daner/daner_mdd_{cohort}.{ancestries}.{build}.{release}.gz"
	log: "logs/sumstats/daner/daner_mdd_{cohort}.{ancestries}.{build}.{release}.log"
	conda: "../envs/meta.yaml" 
	shell: "sh {input.sh} {input.sumstats} {output} {log}"
	
# Convert vcf sumstats to daner
rule vcf2daner:
	input: sumstats="resources/sumstats/vcf_mdd_{cohort}.{ancestries}.{build}.{release}.gz", vcfgwas="resources/vcf/vendor/r-gwasvcf"
	output: "results/sumstats/daner/daner_mdd_{cohort}.{ancestries}.{build}.{release}.gz"
	log: "logs/sumstats/daner/daner_mdd_{cohort}.{ancestries}.{build}.{release}.log"
	conda: "../envs/vcf.yaml"
	script: "../scripts/meta/vcf2daner.R"
	
# for daner files on genome build hg19
rule hg19:
	input: "results/sumstats/daner/daner_mdd_{cohort}.{ancestries}.hg19.{release}.gz"
	output: "results/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{release}.gz"
	log: "logs/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{release}.log"
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
	input: daner="results/sumstats/daner/daner_mdd_{cohort}.{ancestries}.hg38.{release}.gz", chain="resources/liftOver/hg38ToHg19.over.chain"
	output: "results/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{release}.gz"
	log: "logs/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{release}.log"
	conda: "../envs/meta.yaml" 
	script: "../scripts/meta/liftover.R"
	
ruleorder: hg19 > hg38to19

# align to imputation panel
rule align:
	input: daner="results/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.{build}.{release}.gz", ref="results/meta/reference_info"
	output: "results/sumstats/aligned/daner_mdd_{cohort}.{ancestries}.{build}.{release}.aligned.gz"
	log: "logs/sumstats/aligned/daner_mdd_{cohort}.{ancestries}.{build}.{release}.aligned.log"
	conda: "../envs/meta.yaml" 
	script: "../scripts/meta/align.R"

# munge sumstats for ldsc regression
rule meta_ldsc_munge:
	input: sumstats="results/sumstats/aligned/{cohort}.gz", hm3="resources/ldsc/w_hm3.snplist", ldsc=rules.ldsc_install.output
	params:
		prefix="results/ldsc/munged/{cohort}"
	conda: "../envs/ldsc.yaml"
	output: "results/ldsc/munged/{cohort}.sumstats.gz"
	shell: "resources/ldsc/ldsc/munge_sumstats.py --sumstats {input.sumstats} --daner --out {params.prefix} --merge-alleles {input.hm3} --chunksize 500000"
	
# calculate genetic correlation with MDD2
rule meta_ldsc_mdd2:
	input: sumstats="results/ldsc/munged/{cohort}.sumstats.gz", mdd="results/ldsc/munged/daner_mdd_PGC.eur.hg19.wray2018.aligned.sumstats.gz", w_ld=rules.ldsc_unzip_eur_w_ld_chr.output
	params:
		prefix="results/ldsc/rg_mdd/{cohort}"
	conda: "../envs/ldsc.yaml"
	output: "results/ldsc/rg_mdd/{cohort}.log"
	shell: "resources/ldsc/ldsc/ldsc.py --rg {input.sumstats},{input.mdd} --ref-ld-chr {input.w_ld}/ --w-ld-chr {input.w_ld}/ --out {params.prefix}"
	
rg_mdd_logs, = glob_wildcards("results/ldsc/rg_mdd/{cohort}.log")
rule meta_ldsc_mdd2_table:
	input: expand("results/ldsc/rg_mdd/{cohort}.log", cohort=rg_mdd_logs)
	output: "docs/tables/meta_qc_ldsc.txt"
	shell: """tmp=$(mktemp)
	echo -e cohort release rg.mdd2 se gencov > ${{tmp}}.header
	for log in {input}; do 
	sumstats=$(basename $log .log);
	cohort=$(echo $sumstats | awk -F. '{{print $1}}' | awk -F_ '{{print $3}}')
	release=$(echo $sumstats | awk -F. '{{print $4}}')
	gencov=$(cat $log | grep 'Total Observed scale gencov:' | awk '{{print $5}}');
	rg=$(cat $log | grep 'Genetic Correlation:' | awk '{{print $3}}');
	se=$(cat $log | grep 'Genetic Correlation:' | awk '{{print $4}}');
	echo -e $cohort [$release] $rg $se $gencov >> ${{tmp}}.body;
	done;
	cat ${{tmp}}.header > $tmp
	cat ${{tmp}}.body | sort -k 1,2 >> $tmp
	column -t -s' ' $tmp > {output}"""
	
# create reference info file linking to imputation panel
rule refdir:
	output: "results/meta/reference_info"
	log: "logs/meta/reference_info.log"
	shell: "cd results/meta; impute_dirsub --refdir {config[refdir]} --reference_info --outname meta"

# link sumstats files into meta-analysis directory, but also run
# LDSC rg with MDD2
rule meta:
	input: sumstats="results/sumstats/aligned/{cohort}.gz", rg="results/ldsc/rg_mdd/{cohort}.log"
	output: "results/meta/{cohort}.gz"
	log: "logs/meta/{cohort}.log"
	shell: "cp -v {input.sumstats} {output} > {log}"

# Ricopili results dataset list for eur ancestries
rule dataset_eur:
	input: "results/meta/daner_mdd_MDD29.eur.hg19.0120a_rmUKBB.aligned.gz",
	 "results/meta/daner_mdd_23andMe.eur.hg19.v7_2.aligned.gz",
	 "results/meta/daner_mdd_deCODE.eur.hg19.DEPALL_FINAL_WHEAD.aligned.gz",
	 "results/meta/daner_mdd_GenScot.eur.hg19.1215a.aligned.gz",
	 "results/meta/daner_mdd_GERA.eur.hg19.0915a_mds5.aligned.gz",
	 "results/meta/daner_mdd_UKBB.eur.hg19.MD_glm_202010.aligned.gz",
	 "results/meta/daner_mdd_iPSYCH.eur.hg19.2012_HRC.aligned.gz",
	 "results/meta/daner_mdd_iPSYCH.eur.hg19.2015i_HRC.aligned.gz",
	 "results/meta/daner_mdd_FinnGen.eur.hg19.R5_18032020.aligned.gz",
	 "results/meta/daner_mdd_ALSPAC.eur.hg19.12082019.aligned.gz",
	 "results/meta/daner_mdd_Airwave.eur.hg19.0820.aligned.gz",
	 "results/meta/daner_mdd_PBK.eur.hg19.2020.aligned.gz",
	 "results/meta/daner_mdd_ESTBB.eur.hg19.EstBB.aligned.gz",
	 "results/meta/daner_mdd_MoBa.eur.hg19.harvest12.aligned.gz",
	 "results/meta/daner_mdd_MoBa.eur.hg19.harvest24.aligned.gz",
	 "results/meta/daner_mdd_MoBa.eur.hg19.rotterdam1.aligned.gz",
	 "results/meta/daner_mdd_HUNT.eur.hg19.gp_all_20190625.aligned.gz",
	 "results/meta/daner_mdd_HUNT.eur.hg19.hospital_all_20190625.aligned.gz",
	 "results/meta/daner_mdd_STAGE.eur.hg19.MDDdx_fastGWAS.aligned.gz",
	 "results/meta/daner_mdd_PREFECT.eur.hg19.run1.aligned.gz"
	output: "results/meta/dataset_full_eur_v{analysis}"
	log: "logs/meta/dataset_full_eur_v{analysis}.log"
	shell: "for daner in {input}; do echo $(basename $daner) >> {output}; done"
	
# Ricopili results dataset list for eas ancestries
rule dataset_eas:
	input: 
	 "results/meta/daner_mdd_23andMe.eas.hg19.v7_2.aligned.gz",
	 "results/meta/daner_mdd_Taiwan.eas.hg19.20200327.aligned.gz"
	output: "results/meta/dataset_full_eas_v{analysis}"
	log: "logs/meta/dataset_full_eas_v{analysis}.log"
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
analysis_version = ["3.29.15"]
rule postimp_eur:
	input: expand("results/meta/full_eur_v{version}.done", version=analysis_version)
	
rule postimp_eas:
	input: expand("results/meta/full_eas_v{version}.done", version=["3.00.02"])
	
# cohort sets for analysts
cohorts_analyst = ["full", "noUKBB", "no23andMe", "noALSPAC"]

# cohort sets for public
cohorts_public = ["no23andMe"]

