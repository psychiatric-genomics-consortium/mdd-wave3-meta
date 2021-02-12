# QC checks for meta-analysis

# calculate genetic correlation with MDD2 and MDD29
rule meta_ldsc_mdd2:
	input: sumstats="results/sumstats/munged/{cohort}.sumstats.gz", mdd2="results/sumstats/munged/daner_mdd_PGC.eur.hg19.wray2018.aligned.sumstats.gz", mdd29="results/sumstats/munged/daner_mdd_MDD29.eur.hg19.0120a_rmUKBB.aligned.sumstats.gz", w_ld=rules.ldsc_unzip_eur_w_ld_chr.output
	params:
		prefix="results/sumstats/rg_mdd/{cohort}"
	conda: "../envs/ldsc.yaml"
	output: "results/sumstats/rg_mdd/{cohort}.log"
	shell: "resources/ldsc/ldsc/ldsc.py --rg {input.sumstats},{input.mdd2},{input.mdd29} --ref-ld-chr {input.w_ld}/ --w-ld-chr {input.w_ld}/ --out {params.prefix}"
	
rg_mdd_logs, = glob_wildcards("results/sumstats/rg_mdd/{cohort}.log")
rule meta_ldsc_mdd2_table:
	input: expand("results/sumstats/rg_mdd/{cohort}.log", cohort=rg_mdd_logs)
	output: "docs/tables/meta_qc_ldsc.txt"
	shell: """tmp=$(mktemp)
	echo -e cohort release rg.mdd2 se.mdd2 gencov.mdd2 rg.mdd29 se.mdd29 gencov.mdd29 > ${{tmp}}.header
	for log in {input}; do 
	sumstats=$(basename $log .log);
	cohort=$(echo $sumstats | awk -F. '{{print $1}}' | awk -F_ '{{print $3}}');
	release=$(echo $sumstats | awk -F. '{{print $4}}');
	gencovs=$(cat $log | grep 'Total Observed scale gencov:' | awk '{{print $5}}');
	rgs=$(cat $log | grep 'Genetic Correlation:' | awk '{{print $3}}');
	ses=$(cat $log | grep 'Genetic Correlation:' | awk '{{print $4}}');
	gencov_mdd2=$(echo $gencovs | awk '{{print $1}}');
	gencov_mdd29=$(echo $gencovs | awk '{{print $2}}');
	rg_mdd2=$(echo $rgs | awk '{{print $1}}');
	rg_mdd29=$(echo $rgs | awk '{{print $2}}');
	se_mdd2=$(echo $ses | awk '{{print $1}}');
	se_mdd29=$(echo $ses | awk '{{print $2}}');
	echo -e $cohort [$release] $rg_mdd2 $se_mdd2 $gencov_mdd2 $rg_mdd29 $se_mdd29 $gencov_mdd29 >> ${{tmp}}.body;
	done;
	cat ${{tmp}}.header > $tmp
	cat ${{tmp}}.body | sort -k 1,2 >> $tmp
	column -t -s' ' $tmp > {output}"""
	
# LDSC rg between two sets of sumstats
rule meta_ldsc_pairwise_rg:
	input: sumstats1="results/sumstats/munged/daner_mdd_{cohort1}.aligned.sumstats.gz", sumstats2="results/sumstats/munged/daner_mdd_{cohort2}.aligned.sumstats.gz", w_ld=rules.ldsc_unzip_eur_w_ld_chr.output
	params:
		prefix="results/sumstats/rg_pairs/{cohort1},{cohort2}"
	output: "results/sumstats/rg_pairs/{cohort1},{cohort2}.log"
	conda: "../envs/ldsc.yaml"
	shell: "resources/ldsc/ldsc/ldsc.py --rg {input.sumstats1},{input.sumstats2} --ref-ld-chr {input.w_ld}/ --w-ld-chr {input.w_ld}/ --out {params.prefix}"	

# create pairwise lists of all munged sumstats
# first use glob to get a list of all munged sumstats
meta_ldsc_munged_sumstats, = glob_wildcards("results/sumstats/munged/daner_mdd_{cohort}.aligned.sumstats.gz")
# we could use expand() here as the input like:
# 	input: expand("results/ldsc/pairs/{cohort1},{cohort2}.log", cohort1=ldsc_munged_sumstats, cohort2=ldsc_munged_sumstats)
# but this is inefficient because it creates combinations of pairs in all orderings,
# e.g. with "ABC" we get "AA AB AC BA BB BC CA CB CC" but in LDSC --rg A,B and --rg B,A are
# the same and --rg A,A is always 1. Thus we make our own list of desired inputs "AB AC BC" 
# using the sorted combinations sumstats
rule meta_ldsc_sumstats_pairs:
	input: ["results/sumstats/rg_pairs/%s,%s.log" % cohorts for cohorts in list(itertools.combinations(sorted(meta_ldsc_munged_sumstats), r=2))]
	
# merge pairwise LDSC genetic correlation logs into a single table
meta_ldsc_sumstats_pairs_logs, = glob_wildcards("results/sumstats/rg_pairs/{cohorts}.log")

rule meta_ldsc_sumstats_pairs_table:
	input: expand("results/sumstats/rg_pairs/{cohorts}.log", cohorts=meta_ldsc_sumstats_pairs_logs)
	output: "docs/tables/meta_qc_ldsc_pairs.txt"
	conda: "../envs/reports.yaml"
	script: "../scripts/meta/ldsc_pairs_table.R"

# table of alignment checks
align_logs, = glob_wildcards("logs/sumstats/aligned/{cohort}.log")
rule meta_align_qc:
	input: expand("logs/sumstats/aligned/{cohort}.log", cohort=align_logs)
	output: "docs/tables/meta_qc_align.txt"
	conda: "../envs/meta.yaml"
	script: "../scripts/meta/align_qc_table.R"
	
# QC checks
rule meta_qc:
	input: expand("results/meta/dataset_full_eur_v{version}", version=analysis_version),
		"docs/tables/meta_qc_align.txt",
		"docs/tables/meta_qc_ldsc.txt",
		"docs/tables/meta_qc_ldsc_pairs.txt",
		"docs/metaqc.html"	
		
# batch submission example to re-run alignment on LISA:
## $ sbatch -t 4:00:00 --wrap "snakemake -j4 --use-conda --cluster 'sbatch -t 20 -n 16' docs/tables/meta_qc_align.txt --groups hg19=group0 daner=group0 align=group0 meta_align_qc=group1 --group-components group0=6 group1=1"
