# Conduct meta-analysis in Ricopili


#################
#               #
# Meta analyses #
#               #
#################


##
## Cohort lists. List of cohort + subcohort/release pairs
## 
cohorts_eur = [["MDD47", "29w2_18w3_1503"], 
["23andMe", "v7_2_202012"],      
["deCODE", "DEPALL_FINAL_WHEAD"],
["GenScot", "1215a"],            
["GERA", "0915a_mds5"],    
["UKBB", "MD_glm_202012"], 
["iPSYCH", "2012_HRC"],        
["iPSYCH", "2015i_HRC"],       
["FinnGen", "R5_18032020"],      
["ALSPAC", "12082019"],        
["Airwave", "0820"],             
["PBK", "2020"],         
["ESTBB", "EstBB"],          
["MoBa", "harvest12"],     
["MoBa", "harvest24"],     
["MoBa", "rotterdam1"],    
["HUNT", "gp_all_20190625"],
["HUNT", "hospital_all_20190625"],
["STAGE", "MDDdx_saige"],    
["PREFECT", "run1"],             
["AGDS", "202012"],        
["lgic2", "202011"],         
["BASIC", "202011"],        
["BioVU", "Cov_SAIGE_202101"],
["EXCEED", "202010"],          
["MVP", "ICDdep_AllSex_202101"],
["tkda1", "run1"],           
["DBDS", "FINAL202103"],
["SHARE", "godartsshare_842021"]]

cohorts_eas=[["23andMe","v7_2"],
["Taiwan", "20200327"]]

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
	
ruleorder: text2daner > daner

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

# Meta-analysis QC parameters
meta_qc_params = {"maf": 0.01, "info": 0.1, "mac": 20}

	
# create reference info file linking to imputation panel
rule refdir:
	output: "results/meta/reference_info"
	log: "logs/meta/reference_info.log"
	shell: "cd results/meta; impute_dirsub --refdir {config[refdir]} --reference_info --outname meta"

# merged imputation panel SNPs
rule impute_frq2:
	input: ref="results/meta/reference_info"
	output: "results/sumstats/impute_frq2.{ancestries}.rds"
	params:
		maf=meta_qc_params['maf']
	log: "logs/sumstats/impute_frq2.{ancestries}.log"
	conda: "../envs/meta.yaml"
	script: "../scripts/meta/impute_frq2.R"

# align to imputation panel
rule align:
	input: daner="results/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.{build}.{release}.gz", ref="results/sumstats/impute_frq2.{ancestries}.rds", script="scripts/meta/align.R"
	params:
		maf=meta_qc_params['maf'],
		info=meta_qc_params['info'],
		mac=meta_qc_params['mac']
	output: "results/sumstats/aligned/daner_mdd_{cohort}.{ancestries}.{build}.{release}.aligned.gz"
	log: "logs/sumstats/aligned/daner_mdd_{cohort}.{ancestries}.{build}.{release}.aligned.log"
	conda: "../envs/meta.yaml" 
	script: "../scripts/meta/align.R"

# munge sumstats for ldsc regression
rule meta_ldsc_munge:
	input: sumstats="results/sumstats/aligned/{cohort}.gz", hm3="resources/ldsc/w_hm3.snplist", ldsc=rules.ldsc_install.output
	params:
		prefix="results/sumstats/munged/{cohort}"
	conda: "../envs/ldsc.yaml"
	output: "results/sumstats/munged/{cohort}.sumstats.gz"
	shell: "resources/ldsc/ldsc/munge_sumstats.py --sumstats {input.sumstats} --daner --out {params.prefix} --merge-alleles {input.hm3} --chunksize 500000"
	
# extract lists of CPIDs and SNPs from aligned sumstats
rule meta_cpids:
	input: sumstats="results/sumstats/aligned/{cohort}.gz"
	output: "results/sumstats/cpids/{cohort}.cpids.gz"
	log: "logs/sumstats/cpids/{cohort}.log"
	shell: "zcat {input.sumstats} | awk -v OFS='\t' '{{print $1, $2, $3}}' | gzip -c > {output}"

# link sumstats files into meta-analysis directory, but also run
# LDSC rg with MDD2
rule meta:
	input: sumstats="results/sumstats/aligned/{cohort}.gz", rg="results/sumstats/rg_mdd/{cohort}.log"
	output: "results/meta/{cohort}.gz"
	log: "logs/meta/{cohort}.log"
	shell: "cp -v {input.sumstats} {output} > {log}"

# Ricopili results dataset list for eur ancestries
rule dataset_eur:
	input: expand("results/meta/daner_mdd_{cohort}.eur.hg19.{release}.aligned.gz", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur])
	output: "results/meta/dataset_full_eur_v{analysis}"
	log: "logs/meta/dataset_full_eur_v{analysis}.log"
	shell: "for daner in {input}; do echo $(basename $daner) >> {output}; done"
	
# Ricopili results dataset list for eas ancestries
rule dataset_eas:
	input: expand("results/meta/daner_mdd_{cohort}.eas.hg19.{release}.aligned.gz", zip, cohort=[cohort[0] for cohort in cohorts_eas], release=[cohort[1] for cohort in cohorts_eas])
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
analysis_version = ["3.47.24"]
rule postimp_eur:
	input: expand("results/meta/full_eur_v{version}.done", version=analysis_version)
	
rule postimp_eas:
	input: expand("results/meta/full_eas_v{version}.done", version=["3.00.02"])
	
# cohort sets for analysts
cohorts_analyst = ["full", "noUKBB", "no23andMe", "noALSPAC"]

# cohort sets for public
cohorts_public = ["no23andMe"]

# check Ricopili output for complete rows and duplicates
rule postimp_rp:
	input: "results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.gz"
	log: "logs/meta/distribution/{analysis}.rp.log"
	conda: "../envs/meta.yaml"
	output: "results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.rp.gz"
	script: "../scripts/meta/rp.R"
