# Putting together reports to summarize results

# Manhattan plot
rule manhattan_png:
	input: expand("results/distribution/manhattan.nog2.pgc_mdd_full_eur_hg19_v{version}.pdf", version=analysis_version)
	output: "docs/figures/manhattan.nog2.eur.png"
	conda: "../envs/reports.yaml"
	shell: "convert -density 300 {input} -resize 50% {output}"

# Case/control counts
rule cohorts_table:
	input: expand("results/distribution/basic.pgc_mdd_full_eur_hg19_v{version}.num.xls", version=analysis_version)
	output: "docs/tables/cohorts.eur.txt"
	conda: "../envs/reports.yaml"
	script: "../scripts/reports/cohorts_table.py"

# GWAS summary
rule gwas_summary:
	input: manhattan="docs/figures/manhattan.nog2.eur.png", casecontrol=expand("results/distribution/basic.pgc_mdd_full_eur_hg19_v{version}.num.xls", version=analysis_version), snps=expand("results/distribution/daner_pgc_mdd_full_eur_hg19_v{version}.gz.p4.clump.areator.sorted.1mhc.summary", version=analysis_version), cohorts_table="docs/tables/cohorts.eur.txt"
	output: "docs/gwas.md"
	run:
		import pandas as pd
		manhattan_png = os.path.relpath(input.manhattan, "docs")
		#casecontrol = pd.read_excel(input.casecontrol)
		md = """
# MDD3 Meta-analysis

![]({})



""".format(manhattan_png)
		with open(output[0], 'w') as out:
			print(md, file=out)

rule report_metaqc:
	input: meta_qc_align="docs/tables/meta_qc_align.txt", meta_qc_ldsc="docs/tables/meta_qc_ldsc.txt", meta_qc_ldsc_pairs="docs/tables/meta_qc_ldsc_pairs.txt", cohorts_mdd="docs/tables/cohorts_mdd.eur.txt", rmd="scripts/reports/metaqc.Rmd"
	params:
		maf=meta_qc_params['maf'],
		info=meta_qc_params['info']
	output: "docs/metaqc.html"
	conda: "../envs/reports.yaml"
	script: "../scripts/reports/metaqc.Rmd"
	
rule report_metaqc_beta_pca:
	input:  assoc=expand("results/sumstats/metaqc/pruned/full_mdd_{cohort}.eur.hg19.{release}.qc.pruned.assoc", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur]), rmd="scripts/reports/metaqc_pca.Rmd"
	output: "docs/metaqc_pca.html"
	conda: "../envs/reports.yaml"
	script: "../scripts/reports/metaqc_pca.Rmd"
	