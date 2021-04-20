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
	
# prune the sumstats for correlation matrix construction
rule metacarpa_prune:
	input: hm3="resources/ldsc/w_hm3.snplist", bed="resources/1kg/1kg_phase1_all.bed"
	params:
		bed_prefix="resources/1kg/1kg_phase1_all",
		prune_prefix="results/sumstats/metacarpa/1kg.hm3"
	output: "results/sumstats/metacarpa/1kg.hm3.prune.in"
	conda: "../envs/meta.yaml"
	shell: """
	plink \
	--bfile {params.bed_prefix} \
	--extract {input.hm3} \
	--indep-pairwise 1000kb 1 0.1 \
	--maf 0.05 \
	--out {params.prune_prefix}
	"""
	
# extract pruned SNPs from sumstats
rule metacarpa_extract:
	input: assoc="results/sumstats/metacarpa/assoc/mdd_{cohort}.aligned.assoc", snplist="results/sumstats/metacarpa/1kg.hm3.prune.in"
	output: "results/sumstats/metacarpa/pruned/mdd_{cohort}.aligned.pruned.assoc"
	shell: "grep -wFf {input.snplist} {input.assoc} > {output}"
	
rule metacarpa_matrix_eur:
	input: assoc=expand("results/sumstats/metacarpa/pruned/mdd_{cohort}.eur.hg19.{release}.aligned.pruned.assoc", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur]), metacarpa="resources/metacarpa/metacarpa"
	params:
		input_args=lambda wildcards, input: ' '.join(['-I ' + assoc for assoc in input.assoc]),
		output_prefix="results/sumstats/metacarpa/full_eur"
	output: "results/sumstats/metacarpa/full_eur.2300.matrix.txt"
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
	--stop
	"""
	
rule metacarpa_eur:
	input: assoc=expand("results/sumstats/metacarpa/mdd_{cohort}.eur.hg19.{release}.aligned.assoc", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur]), matrix="results/sumstats/metacarpa/full_eur.2300.matrix.txt", metacarpa="resources/metacarpa/metacarpa"
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