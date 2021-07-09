#################
#               #
# Meta analyses #
#               #
#################


##
## Cohort lists. List of cohort + subcohort/release pairs
## 
cohorts_eur = [["MDD49", "29w2_20w3_1504"], 
["23andMe", "v7_2_202012"],      
["deCODE", "DEPALL_FINAL_WHEAD"],
["GenScot", "1215a"],            
["GERA", "0915a_mds5"],    
["UKBB", "MD_glm_202012"], 
["iPSYCH", "2012_HRC"],        
["iPSYCH", "2015i_HRC"],       
["FinnGen", "R5_18032020"],      
["ALSPAC", "27022020"],        
["Airwave", "0820"],             
["PBK", "2020"],         
["ESTBB", "EstBB"],          
["MoBa", "harvest12"],     
["MoBa", "harvest24"],     
["MoBa", "rotterdam1"],
["HUNT", "gp_hospital_metacarpa_20190625"],
["STAGE", "MDDdx_saige"],    
["PREFECT", "run1"],             
["AGDS", "202012"],        
["lgic2", "202011"],         
["BASIC", "202011"],        
["BioVU", "NoCov_SAIGE_051821"],
["EXCEED", "202010"],          
["MVP", "4_0ICDdep_202106"],
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
	shell: "ln -sv {input} {output} > {log}"

# Harmonize names of all summary statistics listed under sumstats in config.yaml
rule sumstats:
	input: expand("resources/sumstats/{sumstats}.gz", sumstats=config["sumstats"])
	
ruleorder: text2daner > daner

# For pre-formatted daner files
rule daner:
	input: "resources/sumstats/daner_{cohort}.gz"
	output: "results/sumstats/daner/daner_{cohort}.gz"
	log: "logs/sumstats/daner/daner_{cohort}.log"
	shell: "ln -sv $(readlink -f {input}) {output} > {log}"
	
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
	shell: "ln -sv $(readlink -f {input}) {output} > {log}"

# download hgIN to hgOUT chain
rule hg_chain:
	input: HTTP.remote("hgdownload.soe.ucsc.edu/goldenPath/hg{from}/liftOver/hg{from}ToHg{to}.over.chain.gz", keep_local=True)
	output: "resources/liftOver/hg{from}ToHg{to}.over.chain"
	run:
		 outputName = os.path.basename(input[0])
		 shell("gunzip -c {input} > {output}")
		 
# download GRCh37/GRCh38 conversion-unstable positions (CUPs) https://github.com/cathaloruaidh/genomeBuildConversionls
rule hg_cups:
	input: HTTP.remote("raw.githubusercontent.com/cathaloruaidh/genomeBuildConversion/master/CUP_FILES/FASTA_BED.ALL_GRCh{build}.novel_CUPs.bed")
	log: "logs/resources/liftOver/GRCh{build}.novel_CUPs.log"
	output: "resources/liftOver/FASTA_BED.ALL_GRCh{build}.novel_CUPs.bed"
	shell: "cp -v {input} {output} > {log}"
	
# liftover hg38 to hg19	
rule hg38to19:
	input: daner="results/sumstats/daner/daner_mdd_{cohort}.{ancestries}.hg38.{release}.gz", chain="resources/liftOver/hg38ToHg19.over.chain", cups="resources/liftOver/FASTA_BED.ALL_GRCh38.novel_CUPs.bed"
	output: "results/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{release}.gz"
	log: "logs/sumstats/hg19/daner_mdd_{cohort}.{ancestries}.hg19.{release}.log"
	conda: "../envs/meta.yaml" 
	script: "../scripts/meta/liftover.R"
	
ruleorder: hg19 > hg38to19

# Meta-analysis QC parameters
meta_qc_params = {"maf": 0.001,
				  "info": 0.1,
				  "mac": 20,
				  "secure_frq": 0.20,
				  "diff_frq": 0.15,
			      "clu_p1": 0.0001,
			      "clu_p2": 0.0001,
				  "clu_r2": 0.1,
			      "clu_kb": 3000,
			      "clu_info": 0.6,
			      "clu_maf": 0.01,
			      "cojo_kb": 50}
	
# create reference info file linking to imputation panel
rule refdir:
	output: "results/meta/reference_info"
	log: "logs/meta/reference_info.log"
	shell: "cd results/meta; impute_dirsub --refdir {config[refdir]} --reference_info --outname meta"

# merged imputation panel SNPs
rule impute_frq2:
	input: ref="results/meta/reference_info", cups="resources/liftOver/FASTA_BED.ALL_GRCh37.novel_CUPs.bed"
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
		secure_frq=meta_qc_params['secure_frq'],
	output: daner="results/sumstats/aligned/daner_mdd_{cohort}.{ancestries}.{build}.{release}.aligned.gz", snp_counts="results/sumstats/aligned/mdd_{cohort}.{ancestries}.{build}.{release}.aligned.txt"
	log: "logs/sumstats/aligned/daner_mdd_{cohort}.{ancestries}.{build}.{release}.aligned.log"
	conda: "../envs/meta.yaml" 
	script: "../scripts/meta/align.R"
	
# apply QC filters
rule filter:
	input: daner="results/sumstats/aligned/daner_mdd_{cohort}.{ancestries}.{build}.{release}.aligned.gz", script="scripts/meta/filter.R"
	params:
		maf=meta_qc_params['maf'],
		info=meta_qc_params['info'],
		mac=meta_qc_params['mac'],
		secure_frq=meta_qc_params['secure_frq'],
		diff_frq=meta_qc_params['diff_frq']
	output: daner="results/sumstats/filtered/daner_mdd_{cohort}.{ancestries}.{build}.{release}.qc.gz", snp_counts="results/sumstats/filtered/mdd_{cohort}.{ancestries}.{build}.{release}.qc.txt"
	log: "logs/sumstats/filtered/daner_mdd_{cohort}.{ancestries}.{build}.{release}.qc.log"
	conda: "../envs/meta.yaml" 
	script: "../scripts/meta/filter.R"
	
# Convert OR to Log-Odds
rule meta_vcf_logOR:
	input: "results/sumstats/filtered/daner_{sumstats}.qc.gz"
	output: "results/sumstats/beta/{sumstats}.gz"
	shell: "gunzip -c {input} | tail -n +2 | awk '{{print $1, $3, $4, $5, log($9), $10, $11, $2, $7, $8}}' | gzip -c > {output}"	
	
# Parsing info for VCF conversion	
rule meta_vcf_daner2vcf_json:
	input: daner="results/sumstats/filtered/daner_mdd_{cohort}.{ancestries}.hg{hg}.{release}.qc.gz"
	output: vcf="results/sumstats/vcf/mdd_{cohort}.{ancestries}.hg{hg}.{release}.json"
	log: "logs/sumstats/vcf/mdd_{cohort}.{ancestries}.hg{hg}.{release}.json.log"
	run:
		with gzip.open(input.daner, 'r') as daner:
			headers = daner.readline().split()
			cohort_cases = headers[5].decode().split('_')[2]
			cohort_controls = headers[6].decode().split('_')[2]
		with open(output.vcf, 'w') as out:
			json.dump({"chr_col": 0,
						"pos_col": 1,
						"ea_col": 2,
						"oa_col": 3,
						"beta_col": 4,
						"se_col": 5,
						"pval_col": 6,
						"snp_col": 7,
						"eaf_col": 8,
						"imp_info_col": 9,
						"delimiter": " ",
						"header": False,
						"build": "{build}".format(build=builds[wildcards.hg]),
						"cohort_cases": cohort_cases,
						"cohort_controls": cohort_controls,
						"id": "mdd.{cohort}.{release}.{ancestries}".format(release=wildcards.release, cohort=wildcards.cohort, ancestries=wildcards.ancestries.upper())},
					 out)

# VCF conversion and harmonisation
rule meta_vcf_daner2vcf:
	input: beta="results/sumstats/beta/mdd_{cohort}.{ancestries}.hg{hg}.{release}.gz", json="results/sumstats/vcf/mdd_{cohort}.{ancestries}.hg{hg}.{release}.json", fasta=lambda wildcards: expand("resources/fasta/human_{build}.fasta", build=builds[wildcards.hg].lower()), fai=lambda wildcards: expand("resources/fasta/human_{build}.fasta.fai", build=builds[wildcards.hg].lower()), gwas2vcf=rules.vcf_install_gwas2vcf.output
	output: vcf="results/sumstats/vcf/mdd_{cohort}.{ancestries}.hg{hg}.{release}.vcf.gz", tbi="results/sumstats/vcf/mdd_{cohort}.{ancestries}.hg{hg}.{release}.vcf.gz.tbi"
	log: "logs/sumstats/vcf/mdd_{cohort}.{ancestries}.hg{hg}.{release}.log"
	conda: "../envs/vcf.yaml"
	shell: "python resources/vendor/gwas2vcf/main.py --out {output.vcf} --data {input.beta} --json {input.json} --ref {input.fasta} > {log}"

# Merge VCF files
rule meta_vcf_merge_eur:
	input: expand("results/sumstats/vcf/mdd_{cohort}.eur.hg19.{release}.vcf.gz", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur])
	conda: "../envs/vcf.yaml"
	output: "results/sumstats/mdd_cohorts_eur.vcf.gz"
	shell: "bcftools merge -O z -o {output} {input}"
	
rule meta_vcf_merge_eas:
	input: expand("results/sumstats/vcf/mdd_{cohort}.eas.hg19.{release}.vcf.gz", zip, cohort=[cohort[0] for cohort in cohorts_eas], release=[cohort[1] for cohort in cohorts_eas])
	conda: "../envs/vcf.yaml"
	output: "results/sumstats/mdd_cohorts_eas.vcf.gz"
	shell: "bcftools merge -O z -o {output} {input}"
	

# munge sumstats for ldsc regression
rule meta_ldsc_munge:
	input: sumstats="results/sumstats/filtered/{cohort}.gz", hm3="resources/ldsc/w_hm3.snplist", ldsc=rules.ldsc_install.output
	params:
		prefix="results/sumstats/munged/{cohort}"
	conda: "../envs/ldsc.yaml"
	output: "results/sumstats/munged/{cohort}.sumstats.gz"
	shell: "resources/ldsc/ldsc/munge_sumstats.py --sumstats {input.sumstats} --daner --out {params.prefix} --merge-alleles {input.hm3} --chunksize 500000"
	
# calculate observed scale h2
rule meta_ldsc_h2:
	input: sumstats="results/sumstats/munged/{cohort}.sumstats.gz", w_ld=rules.ldsc_unzip_eur_w_ld_chr.output
	params:
		prefix="results/sumstats/h2/{cohort}"
	conda: "../envs/ldsc.yaml"
	output: "results/sumstats/h2/{cohort}.log"
	shell: "resources/ldsc/ldsc/ldsc.py --h2 {input.sumstats} --ref-ld-chr {input.w_ld}/ --w-ld-chr {input.w_ld}/ --out {params.prefix}"
	
rule meta_ldsc_h2_table:
	input: expand("results/sumstats/h2/daner_mdd_{cohort}.eur.hg19.{release}.qc.log", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur])
	output: "docs/tables/meta_ldsc_h2.txt"
	shell: """
tmp=$(mktemp)
echo -e "cohort\trelease\th2_obs\th2_obs_se\tlambdaGC\tmeanChisq\tintercept\tintercept_se" > ${{tmp}}.header
for log in {input}; do
sumstats=$(basename $log .log);
cohort=$(echo $sumstats | awk -F. '{{print $1}}' | awk -F_ '{{print $3}}');
release=$(echo $sumstats | awk -F. '{{print $4}}');
h2_entry=$(cat $log | grep 'Total Observed scale h2:' || true)
lambda_entry=$(cat $log | grep 'Lambda GC:' || true)
chisq_entry=$(cat $log | grep 'Mean Chi^2:' || true)
intercept_entry=$(cat $log | grep 'Intercept:' || true)
ratio_entry=$(cat $log | grep 'Ratio:' || true)
h2=$(echo $h2_entry | awk '{{print $5}}');
h2_se=$(echo $h2_entry | awk '{{print $6}}' | sed -e 's/[()]//g');
lambdagc=$(echo $lambda_entry | awk '{{print $3}}');
chisq=$(echo $chisq_entry | awk '{{print $3}}');
intercept=$(echo $intercept_entry | awk '{{print $2}}');
intercept_se=$(echo $intercept_entry | awk '{{print $3}}' | sed -e 's/[()]//g');
echo -e "$cohort\t[$release]\t$h2\t$h2_se\t$lambdagc\t$chisq\t$intercept\t$intercept_se"  >> ${{tmp}}.body;
done;
cat ${{tmp}}.header > {output};
cat ${{tmp}}.body | sort -k 1,2 >> {output}"""

# extract lists of CPIDs and SNPs from aligned sumstats
rule meta_cpids:
	input: sumstats="results/sumstats/aligned/{cohort}.gz"
	output: "results/sumstats/cpids/{cohort}.cpids.gz"
	log: "logs/sumstats/cpids/{cohort}.log"
	shell: "zcat {input.sumstats} | awk -v OFS='\t' '{{print $1, $2, $3}}' | gzip -c > {output}"

