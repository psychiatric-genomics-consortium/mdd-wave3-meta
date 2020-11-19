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
