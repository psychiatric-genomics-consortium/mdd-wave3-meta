# Meta-analysis using METACARPA https://github.com/hmgu-itg/metacarpa/

	
# list of cohorts/releases
# split into groups of 16 or fewer (==bitmask size in metacarpa)
cohorts_eur_16 = {
'group1': [['23andMe', 'v7_2_202012'],
['MVP', 'ICDdep_AllSex_202101'],
['UKBB', 'MD_glm_202012'],
['MDD49', '29w2_20w3_1504'],
['AGDS', '202012'],
['BioVU', 'Cov_SAIGE_202101'],
['GERA', '0915a_mds5'],
['PBK', '2020'],
['Airwave', '0820'],
['SHARE', 'godartsshare_842021'],
['GenScot', '1215a'],
['EXCEED', '202010'],
['ALSPAC', '12082019']],
'group2': [['FinnGen', 'R5_18032020'],
['ESTBB', 'EstBB'],
['deCODE', 'DEPALL_FINAL_WHEAD'],
['iPSYCH', '2012_HRC'],
['DBDS', 'FINAL202103'],
['iPSYCH', '2015i_HRC'],
['HUNT', 'gp_all_20190625'],
['HUNT', 'hospital_all_20190625'],
['PREFECT', 'run1'],
['BASIC', '202011'],
['lgic2', '202011'],
['tkda1', 'run1'],
['MoBa', 'harvest12'],
['MoBa', 'rotterdam1'],
['STAGE', 'MDDdx_saige'],
['MoBa', 'harvest24']]}

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

# SNPs common to cohort groups	
rule metacarpa_snps_eur:
	input: lambda wildcards: expand("results/sumstats/metacarpa/snps/mdd_{cohort}.eur.hg19.{release}.aligned.snplist", zip, cohort=[cohort[0] for cohort in cohorts_eur_16[wildcards.group]], release=[cohort[1] for cohort in cohorts_eur_16[wildcards.group]])
	params:
		n_cohorts=len(cohorts_eur),
		mask_min=30000
	conda: "../envs/meta.yaml"
	output: "results/sumstats/metacarpa/{group}_eur.snplist"
	script: "../scripts/meta/metacarpa_snplist.py"

# extract common SNPs	
rule metacarpa_extract:
	input: assoc="results/sumstats/metacarpa/assoc/mdd_{cohort}.{ancestries}.hg19.{release}.aligned.assoc", snplist="results/sumstats/metacarpa/{group}_eur.snplist"
	output: "results/sumstats/metacarpa/extract/{group}_mdd_{cohort}.{ancestries}.hg19.{release}.aligned.assoc"
	shell: "grep -wFf {input.snplist} {input.assoc} > {output}"

# prune the sumstats for correlation matrix construction
rule metacarpa_prune:
	input: hm3="resources/ldsc/w_hm3.snplist", snplist="results/sumstats/metacarpa/{group}_{ancestries}.snplist", bed="resources/1kg/1kg_phase1_all.bed"
	params:
		bed_prefix="resources/1kg/1kg_phase1_all",
		prune_prefix="results/sumstats/metacarpa/{group}.1kg.{ancestries}"
	output: "results/sumstats/metacarpa/{group}.1kg.{ancestries}.prune.in"
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
	input: assoc="results/sumstats/metacarpa/extract/{group}_mdd_{cohort}.{ancestries}.hg19.{release}.aligned.assoc", snplist="results/sumstats/metacarpa/{group}.1kg.{ancestries}.prune.in"
	output: "results/sumstats/metacarpa/pruned/{group}_mdd_{cohort}.{ancestries}.hg19.{release}.aligned.pruned.assoc"
	shell: "grep -wFf {input.snplist} {input.assoc} > {output}"


	
rule metacarpa_matrix_eur:
	input: assoc=lambda wildcards: expand("results/sumstats/metacarpa/pruned/{{group}}_mdd_{cohort}.eur.hg19.{release}.aligned.pruned.assoc", zip, cohort=[cohort[0] for cohort in cohorts_eur_16[wildcards.group]], release=[cohort[1] for cohort in cohorts_eur_16[wildcards.group]]), metacarpa="resources/metacarpa/metacarpa"
	params:
		input_args=lambda wildcards, input: ' '.join(['-I ' + assoc for assoc in input.assoc]),
		output_prefix="results/sumstats/metacarpa/{group}_eur"
	output: "results/sumstats/metacarpa/{group}_eur.matrix.txt"
	shell: """
	{input.metacarpa} \
	{params.input_args} \
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
	input: assoc=lambda wildcards: expand("results/sumstats/metacarpa/extract/{{group}}_mdd_{cohort}.eur.hg19.{release}.aligned.assoc", zip, cohort=[cohort[0] for cohort in cohorts_eur_16[wildcards.group]], release=[cohort[1] for cohort in cohorts_eur_16[wildcards.group]]), matrix="results/sumstats/metacarpa/{group}_eur.matrix.txt", metacarpa="resources/metacarpa/metacarpa"
	params:
		input_args=lambda wildcards, input: ' '.join(['-I ' + assoc for assoc in input.assoc])
	output: "results/meta/metacarpa/pgc_mdd_{group}_eur_hg19_v{version}.mc.txt"
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
	--size-col 10
	"""

# Run meta-carpa for both groups	
rule metacarpa_eur_analyze:
	input: expand("results/meta/metacarpa/pgc_mdd_{group}_eur_hg19_v{version}.mc.txt", version=analysis_version, group=['group1', 'group2'])
	
# format into assoc files for input back into metacarpa
#  1	rsid
#  2	chr:pos
#  3	effect_allele
#  4	neffect_allele
#  5	effect_allele_frequency
#  6	effects
#  7	beta
#  8	se
#  9	z
# 10	z_se
# 11	p_wald
# 12	p_corrected
# 13    p_stouffer
# 14	n
rule metacarpa_group_assoc:
	input: "results/meta/metacarpa/pgc_mdd_{group}_eur_hg19_v{version}.mc.txt"
	output: "results/meta/metacarpa/pgc_mdd_{group}_eur_hg19_v{version}.mc.assoc"
	shell: """
	cat {input} | awk 'NR > 1 {{sub(":", " ", $2); print $2, $3, $4, $11, $7, $8, $5, $1, $14}}' > {output}
	"""
	
rule metacarpa_eur_meta:
	input: group1="results/meta/metacarpa/pgc_mdd_group1_eur_hg19_v{version}.mc.assoc", group2="results/meta/metacarpa/pgc_mdd_group2_eur_hg19_v{version}.mc.assoc", metacarpa="resources/metacarpa/metacarpa"
	output: "results/meta/metacarpa/pgc_mdd_full_eur_hg19_v{version}.mc.meta.txt"
	shell: """
	{input.metacarpa} -I {input.group1} -I {input.group2} \
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
	--size-col 10
	"""
	
rule metacarpa_group_eur_analyze:
	input: expand("results/meta/metacarpa/pgc_mdd_full_eur_hg19_v{version}.mc.meta.txt", version=analysis_version)
	
rule metacarpa_group_daner:
	input: "results/meta/metacarpa/pgc_mdd_full_{analysis}.mc.meta.txt"
	output: "results/meta/metacarpa/daner_pgc_mdd_full_{analysis}.mc.meta.gz"
	conda: "../envs/meta.yaml"
	script: "../scripts/meta/metacarpa_daner.R"
