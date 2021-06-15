# Conduct meta-analysis in Ricopili

# link sumstats files into meta-analysis directory, but also run
# LDSC rg with MDD2
rule meta:
	input: sumstats="results/sumstats/aligned/{cohort}.gz", rg="results/sumstats/rg_mdd/{cohort}.log"
	output: "results/meta/{cohort}.gz"
	log: "logs/meta/{cohort}.log"
	shell: "ln -sv $(readlink -f {input.sumstats}) {output} > {log}"

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
	
# inputs for postimp_rp
rule postimp_rp_all:
	input: expand("results/meta/distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.rp.gz", cohorts=cohorts_analyst, ancestries=['eur'], version=analysis_version)