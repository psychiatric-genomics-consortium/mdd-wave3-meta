# Meta-analysis using METACARPA https://github.com/hmgu-itg/metacarpa/

# Download the METACARPA binary
rule metacarpa_download:
	input: HTTP.remote("https://github.com/hmgu-itg/metacarpa/releases/download/1.0.1/metacarpa", keep_local=False)
	output: "resources/metacarpa/metacarpa"
	shell: "cp {input} {output}"
	
# format aligned sumstats to METACARPA assoc input
# CHR POS A1 A2 P BETA SE AF Neff
# Neff = (4 * Ncases * Ncontrols) / (Ncases + Ncontrols)
rule metacarpa_sumstats:
	input: "results/sumstats/aligned/daner_mdd_{cohort}.aligned.gz"
	output: "results/sumstats/metacarpa/assoc/mdd_{cohort}.aligned.assoc"
	shell: "zcat {input} | awk -v neff=$(zcat {input} | head -n 1 | awk '{{print $6, $7}}' | sed 's/_/ /g' | awk '{{print (4*$3*$6)/($3+$6)}}') 'NR > 1 {{print $1, $3, $4, $5, $11, log($9), $10, $7, $2, neff}}' > {output}"
	
# get list of SNPs in each cohort
rule metacarpa_snps:
	input: "results/sumstats/metacarpa/assoc/mdd_{cohort}.aligned.assoc"
	output: "results/sumstats/metacarpa/snps/mdd_{cohort}.aligned.snplist"
	shell: "cat {input} | awk '{{print $9}}' | sort > {output}"

# SNPs common to all cohorts	
rule metacarpa_snps_eur:
	input: expand("results/sumstats/metacarpa/snps/mdd_{cohort}.eur.hg19.{release}.aligned.snplist", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur])
	params:
		n_cohorts=len(cohorts_eur)
	output: "results/sumstats/metacarpa/full_eur.snplist"
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

# extract common SNPs	
rule metacarpa_extract:
	input: assoc="results/sumstats/metacarpa/assoc/mdd_{cohort}.{ancestries}.hg19.{release}.aligned.assoc", snplist="results/sumstats/metacarpa/full_eur.snplist"
	output: "results/sumstats/metacarpa/extract/mdd_{cohort}.{ancestries}.hg19.{release}.aligned.assoc"
	shell: "grep -wFf {input.snplist} {input.assoc} > {output}"

# prune the sumstats for correlation matrix construction
rule metacarpa_prune:
	input: hm3="resources/ldsc/w_hm3.snplist", snplist="results/sumstats/metacarpa/full_{ancestries}.snplist", bed="resources/1kg/1kg_phase1_all.bed"
	params:
		bed_prefix="resources/1kg/1kg_phase1_all",
		prune_prefix="results/sumstats/metacarpa/1kg.{ancestries}"
	output: "results/sumstats/metacarpa/1kg.{ancestries}.prune.in"
	conda: "../envs/meta.yaml"
	shell: """
	plink \
	--bfile {params.bed_prefix} \
	--extract {input.snplist} \
	--indep-pairwise 500kb 1 0.2 \
	--maf 0.05 \
	--out {params.prune_prefix}
	"""
	
# extract pruned SNPs from sumstats
rule metacarpa_prune_extract:
	input: assoc="results/sumstats/metacarpa/extract/mdd_{cohort}.{ancestries}.hg19.{release}.aligned.assoc", snplist="results/sumstats/metacarpa/1kg.{ancestries}.prune.in"
	output: "results/sumstats/metacarpa/pruned/mdd_{cohort}.{ancestries}.hg19.{release}.aligned.pruned.assoc"
	shell: "grep -wFf {input.snplist} {input.assoc} > {output}"
	
# sample-size sorted list of cohorts
cohorts_eur_n = [['23andMe', 'v7_2_202012'],
['MVP', 'ICDdep_AllSex_202101'],
['UKBB', 'MD_glm_202012'],
['ESTBB', 'EstBB'],
['MDD49', '29w2_20w3_1504'],
['FinnGen', 'R5_18032020'],
['deCODE', 'DEPALL_FINAL_WHEAD'],
['iPSYCH', '2012_HRC'],
['DBDS', 'FINAL202103'],
['AGDS', '202012'],
['iPSYCH', '2015i_HRC'],
['BioVU', 'Cov_SAIGE_202101'],
['GERA', '0915a_mds5'],
['HUNT', 'gp_all_20190625'],
['PBK', '2020']]

# ['HUNT', 'hospital_all_20190625']],
# ['Airwave', '0820'],
# ['PREFECT', 'run1'],
# ['SHARE', 'godartsshare_842021'],
# ['BASIC', '202011'],
# ['GenScot', '1215a'],
# ['lgic2', '202011'],
# ['tkda1', 'run1'],
# ['MoBa', 'harvest12'],
# ['EXCEED', '202010'],
# ['MoBa', 'rotterdam1'],
# ['ALSPAC', '12082019'],
# ['STAGE', 'MDDdx_saige'],
# ['MoBa', 'harvest24']]

	
rule metacarpa_matrix_eur:
	input: assoc=expand("results/sumstats/metacarpa/pruned/mdd_{cohort}.eur.hg19.{release}.aligned.pruned.assoc", zip, cohort=[cohort[0] for cohort in cohorts_eur_n], release=[cohort[1] for cohort in cohorts_eur_n]), metacarpa="resources/metacarpa/metacarpa"
	params:
		input_args=lambda wildcards, input: ' '.join(['-I ' + assoc for assoc in input.assoc]),
		output_prefix="results/sumstats/metacarpa/full_eur"
	output: "results/sumstats/metacarpa/full_eur.matrix.txt"
	shell: """
	{input.metacarpa} {params.input_args} \
	--output {params.output_prefix} \
	--sep ' ' \
	--chr-col 1 \
	--pos-col 2 \
	--a1-col 3 \
	--a2-col 4 \
	--pval-col 5 \
	--beta-col 6 \
	--se-col 7 \
	--af-col 8 \
	--id-col 9 \
	--size-col 10 \
	--stop;
	mv {params.output_prefix}.*.matrix.txt {output}
	"""
	
rule metacarpa_eur:
	input: assoc=expand("results/sumstats/metacarpa/extract/mdd_{cohort}.eur.hg19.{release}.aligned.assoc", zip, cohort=[cohort[0] for cohort in cohorts_eur_n], release=[cohort[1] for cohort in cohorts_eur_n]), matrix="results/sumstats/metacarpa/full_eur.matrix.txt", metacarpa="resources/metacarpa/metacarpa"
	params:
		input_args=lambda wildcards, input: ' '.join(['-I ' + assoc for assoc in input.assoc])
	output: "results/meta/metacarpa/pgc_mdd_full_eur_hg19_v{version}.mc.txt"
	shell: """
	{input.metacarpa} {params.input_args} \
	--output {output} \
	--sep ' ' \
	--chr-col 1 \
	--pos-col 2 \
	--a1-col 3 \
	--a2-col 4 \
	--pval-col 5 \
	--beta-col 6 \
	--se-col 7 \
	--af-col 8 \
	--id-col 9 \
	--size-col 10 \
	--matrix {input.matrix}
	"""
	
rule metacarpa_eur_analyze:
	input: expand("results/meta/metacarpa/pgc_mdd_full_eur_hg19_v{version}.mc.txt", version=analysis_version)