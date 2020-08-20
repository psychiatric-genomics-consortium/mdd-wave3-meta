# Putting together reports to summarize results

# Manhattan plot
rule manhattan_png:
	input: expand("results/distribution/manhattan.nog2.pgc_mdd_full_eur_hg19_v{version}.pdf", version=analysis_version)
	output: "docs/figures/manhattan.nog2.eur.png"
	conda: "../envs/reports.yaml"
	shell: "convert -density 300 {input} -resize 50% {output}"