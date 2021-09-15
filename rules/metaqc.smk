# QC checks for meta-analysis

# calculate genetic correlation with MDD2 and MDD29
rule meta_ldsc_mdd2:
	input: sumstats="results/sumstats/munged/{cohort}.sumstats.gz", mdd2="results/sumstats/munged/daner_mdd_PGC.eur.hg19.wray2018.qc.sumstats.gz", mdd29="results/sumstats/munged/daner_mdd_MDD29.eur.hg19.0120a_rmUKBB.qc.sumstats.gz", w_ld=rules.ldsc_unzip_eur_w_ld_chr.output
	params:
		prefix="results/sumstats/rg_mdd/{cohort}"
	conda: "../envs/ldsc.yaml"
	output: "results/sumstats/rg_mdd/{cohort}.log"
	shell: "resources/ldsc/ldsc/ldsc.py --rg {input.sumstats},{input.mdd2},{input.mdd29} --ref-ld-chr {input.w_ld}/ --w-ld-chr {input.w_ld}/ --out {params.prefix}"
	
rg_mdd_logs, = glob_wildcards("results/sumstats/rg_mdd/{cohort}.log")
rule meta_ldsc_mdd2_table:
	input: expand("results/sumstats/rg_mdd/{cohort}.log", cohort=rg_mdd_logs)
	output: "docs/tables/meta_qc_ldsc.txt"
	shell: """
tmp=$(mktemp)
echo -e "cohort\trelease\th2_obs\th2_obs_se\trg.mdd2\tse.mdd2\tgencov.mdd2\tgencov_se.mdd2\tz1z2.mdd2\trg.mdd29\tse.mdd29\tgencov.mdd29\tgencov_se.mdd29\tz1z2.mdd29" > ${{tmp}}.header
for log in {input}; do
sumstats=$(basename $log .log);
cohort=$(echo $sumstats | awk -F. '{{print $1}}' | awk -F_ '{{print $3}}');
release=$(echo $sumstats | awk -F. '{{print $4}}');
h2_entry=$(cat $log | awk '/Heritability of phenotype 1/,/Heritability of phenotype 2/' | grep 'Total Observed scale h2:' || true)
gencov_mdd2_entry=$(cat $log | awk '/rg for phenotype 2/,/rg for phenotype 3/' | grep 'Total Observed scale gencov:' || true);
z1z2_mdd2_entry=$(cat $log | awk '/rg for phenotype 2/,/rg for phenotype 3/' | grep 'Mean z1' || true);
rg_mdd2_entry=$(cat $log | awk '/rg for phenotype 2/,/rg for phenotype 3/' | grep 'Genetic Correlation:' || true);
gencov_mdd29_entry=$(cat $log | awk '/rg for phenotype 3/,/Summary of Genetic/' | grep 'Total Observed scale gencov:' || true);
z1z2_mdd29_entry=$(cat $log | awk '/rg for phenotype 3/,/Summary of Genetic/' | grep 'Mean z1' || true);
rg_mdd29_entry=$(cat $log | awk '/rg for phenotype 3/,/Summary of Genetic/' | grep 'Genetic Correlation:' || true);
h2=$(echo $h2_entry | awk '{{print $5}}');
h2_se=$(echo $h2_entry | awk '{{print $6}}' | sed -e 's/[()]//g');
gencov_mdd2=$(echo $gencov_mdd2_entry | awk '{{print $5}}');
gencov_se_mdd2=$(echo $gencov_mdd2_entry | awk '{{print $6}}'| sed -e 's/[()]//g');
z1z2_mdd2=$(echo $z1z2_mdd2_entry | awk '{{print $3}}');
rg_mdd2=$(echo $rg_mdd2_entry | awk '{{print $3}}');
se_mdd2=$(echo $rg_mdd2_entry | awk '{{print $4}}'| sed -e 's/[()]//g');
gencov_mdd29=$(echo $gencov_mdd29_entry | awk '{{print $5}}');
gencov_se_mdd29=$(echo $gencov_mdd29_entry | awk '{{print $6}}'| sed -e 's/[()]//g');
z1z2_mdd29=$(echo $z1z2_mdd29_entry | awk '{{print $3}}');
rg_mdd29=$(echo $rg_mdd29_entry | awk '{{print $3}}');
se_mdd29=$(echo $rg_mdd29_entry | awk '{{print $4}}'| sed -e 's/[()]//g');
echo -e "$cohort\t[$release]\t$h2\t$h2_se\t$rg_mdd2\t$se_mdd2\t$gencov_mdd2\t$gencov_se_mdd2\t$z1z2_mdd2\t$rg_mdd29\t$se_mdd29\t$gencov_mdd29\t$gencov_se_mdd29\t$z1z2_mdd29"  >> ${{tmp}}.body;
done;
cat ${{tmp}}.header > {output};
cat ${{tmp}}.body | sort -k 1,2 >> {output}"""
	
# LDSC rg between two sets of sumstats
rule meta_ldsc_pairwise_rg:
	input: sumstats1="results/sumstats/munged/daner_mdd_{cohort1}.qc.sumstats.gz", sumstats2="results/sumstats/munged/daner_mdd_{cohort2}.qc.sumstats.gz", w_ld=rules.ldsc_unzip_eur_w_ld_chr.output
	params:
		prefix="results/sumstats/rg_pairs/{cohort1},{cohort2}"
	output: "results/sumstats/rg_pairs/{cohort1},{cohort2}.log"
	conda: "../envs/ldsc.yaml"
	shell: "resources/ldsc/ldsc/ldsc.py --rg {input.sumstats1},{input.sumstats2} --ref-ld-chr {input.w_ld}/ --w-ld-chr {input.w_ld}/ --out {params.prefix}"	

# create pairwise lists of all munged sumstats
# first use glob to get a list of all munged sumstats
meta_ldsc_munged_sumstats, = glob_wildcards("results/sumstats/munged/daner_mdd_{cohort}.qc.sumstats.gz")
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
	input: aligned=expand("results/sumstats/aligned/mdd_{cohort}.{ancestry}.hg19.{release}.aligned.txt", zip, cohort=[cohort[0] for cohort in cohorts_eur + cohorts_eas], ancestry=['eur'] * len(cohorts_eur) + ['eas'] * len(cohorts_eas), release=[cohort[1] for cohort in cohorts_eur + cohorts_eas]), filtered=expand("results/sumstats/filtered/mdd_{cohort}.{ancestry}.hg19.{release}.qc.txt", zip, cohort=[cohort[0] for cohort in cohorts_eur + cohorts_eas], ancestry=['eur'] * len(cohorts_eur) + ['eas'] * len(cohorts_eas), release=[cohort[1] for cohort in cohorts_eur + cohorts_eas])
	output: "docs/tables/meta_qc_align.txt"
	conda: "../envs/meta.yaml"
	script: "../scripts/meta/align_qc_table.R"
	
# get list of SNPs in each cohort
rule meta_snps:
	input: "results/sumstats/filtered/daner_mdd_{cohort}.qc.gz"
	output: "results/sumstats/snps/mdd_{cohort}.qc.snplist"
	shell: "zcat {input} | awk 'NR > 1 {{print $2}}' | sort > {output}"
	
# SNPs common to all cohorts	
rule meta_snps_eur:
	input: expand("results/sumstats/snps/mdd_{cohort}.eur.hg19.{release}.qc.snplist", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur])
	params:
		n_cohorts=len(cohorts_eur)
	output: "results/sumstats/full_eur.snplist"
	run:
		# get SNPs from first file and make a set
		with open(input[0]) as f:
			snps = set(f.read().splitlines())
		# open all SNP lists one by one and intersect them with the snps set
		for snplist in input:
			print(snplist)
			with open(snplist) as f:
				snps = set.intersection(snps, set(f.read().splitlines()))
			print(len(snps))
		# output intersection of all SNPs
		with open(output[0], 'w') as out:
			out.writelines(map(lambda snp: snp+'\n', list(snps))) 
			
rule meta_qc_prune:
	input: snplist="results/sumstats/{group}_{ancestries}.snplist", bed="resources/1kg/1kg_phase1_all.bed"
	params:
		bed_prefix="resources/1kg/1kg_phase1_all",
		prune_prefix="results/sumstats/{group}.1kg.{ancestries}"
	output: "results/sumstats/{group}.1kg.{ancestries}.prune.in"
	conda: "../envs/meta.yaml"
	shell: """
	plink \
	--bfile {params.bed_prefix} \
	--extract {input.snplist} \
	--indep-pairwise 500kb 1 0.2 \
	--maf 0.05 \
	--out {params.prune_prefix} \
	--memory 2000
	"""
	
		
# extract pruned SNPs from sumstats
rule meta_qc_prune_extract:
	input: assoc= "results/sumstats/filtered/daner_mdd_{cohort}.{ancestries}.hg19.{release}.qc.gz", snplist="results/sumstats/{group}.1kg.{ancestries}.prune.in"
	output: "results/sumstats/metaqc/pruned/{group}_mdd_{cohort}.{ancestries}.hg19.{release}.qc.pruned.assoc"
	shell: "zcat {input.assoc} | awk '{{print $1, $2, $3, $7, $9, $10}}' | grep -wFf {input.snplist} > {output}"
	

	
# extract top SNPs from sumstats
rule meta_qc_snps_clumped:
	input: expand("results/distribution/daner_pgc_mdd_{{group}}_{{ancestries}}_hg19_v{version}.gz.p4.clump.areator.sorted.1mhc", version=analysis_version)
	output: "results/sumstats/metaqc/{group}.{ancestries}.clumped.snplist"
	shell: "cat {input} | awk 'NR > 1 {{print $1}}' > {output}"

rule meta_qc_clumped_extract:
	input: assoc= "results/sumstats/filtered/daner_mdd_{cohort}.{ancestries}.hg19.{release}.qc.gz", snplist="results/sumstats/metaqc/{group}.{ancestries}.clumped.snplist"
	output: "results/sumstats/metaqc/clumped/{group}_mdd_{cohort}.{ancestries}.hg19.{release}.qc.clumped.assoc"
	shell: "zcat {input.assoc} | awk '{{print $1, $2, $3, $7, $9, $10}}' | grep -wFf {input.snplist} > {output}"
	
# QC checks
rule meta_qc:
	input: expand("results/meta/dataset_full_eur_v{version}", version=analysis_version),
		"docs/tables/meta_qc_align.txt",
		"docs/tables/meta_qc_ldsc.txt",
		"docs/tables/meta_qc_ldsc_pairs.txt",
		"docs/metaqc.html"	
		
# batch submission example to re-run alignment on LISA:
## $ sbatch -t 4:00:00 --wrap "snakemake -j4 --use-conda --cluster 'sbatch -t 20 -n 16' docs/tables/meta_qc_align.txt --groups text2daner=group0 hg19=group0 daner=group0 align=group0 filter=group0 meta_align_qc=group1 --group-components group0=4 group1=1"
