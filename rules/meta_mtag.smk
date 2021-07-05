# Meta-analysis using MTAG 

# Install mtag
rule mtag_install:
    output: repo=directory("resources/mtag")
	conda: "../envs/mtag.yaml"
	shell: "git clone https://github.com/omeed-maghzian/mtag.git {output}"
	
# format summary statistics
# snpid, chr, bpos, a1, a2, freq, beta, se, pval, n
# case/control Ns are pulled from the header and converted to effective N
rule mtag_sumstats:
	input: "results/sumstats/filtered/daner_mdd_{cohort}.qc.gz"
	output: "results/sumstats/mtag/mdd_{cohort}.qc.txt"
	shell: """
	zcat {input} | awk -v Neff=$(zcat {input} | head -n 1 | awk '{{print $6, $7}}' | sed "s/_/ /g" | awk '{{print (4*$3*$6)/($3 + $6)}}') '{{if(NR == 1) {{print "snpid", "chr", "bpos", "a1", "a2", "freq", "beta", "se", "pval", "n"}} else {{if(NF == 19) {{print $2, $1, $3, $4, $5, $7, log($9), $10, $11, 2*$19}} else {{print $2, $1, $3, $4, $5, $7, log($9), $10, $11, Neff}}}}}}' > {output}
	"""

# run mtag across sumstats
rule mtag_eur:
	input: sumstats=expand("results/sumstats/mtag/mdd_{cohort}.eur.hg19.{release}.qc.txt", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur]), mtag=rules.mtag_install.output
	params:
		input_args=lambda wildcards, input: ','.join(input.sumstats),
		prefix="results/meta/mtag/mdd_eur"
	output: "results/meta/mtag/mdd_eur_trait_1.txt"
	conda: "../envs/mtag.yaml"
	shell: "python resources/mtag/mtag.py --sumstats {params.input_args} --use_beta_se --beta_name beta --se_name se --force --out {params.prefix} --stream_stdout"
	
# run mtag meta-analysis (assume all are measures of the same trait)
rule mtag_meta_eur:
	input: sumstats=expand("results/sumstats/mtag/mdd_{cohort}.eur.hg19.{release}.qc.txt", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur]), mtag=rules.mtag_install.output
	params:
		input_args=lambda wildcards, input: ','.join(input.sumstats),
		prefix="results/meta/mtag/mdd_eur_meta"
	output: "results/meta/mtag/mdd_eur_meta.txt"
	conda: "../envs/mtag.yaml"
	shell: "python resources/mtag/mtag.py --sumstats {params.input_args} --use_beta_se --beta_name beta --se_name se --force --out {params.prefix} --stream_stdout --perfect_gencov --equal_h2 --meta_format"
