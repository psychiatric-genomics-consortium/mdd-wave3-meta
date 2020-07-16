configfile: "config.yaml"

rule daner_link:
	input: lambda wildcards: config["daner"][wildcards.cohort]
	output: "daner/{cohort}.gz"
	shell: "ln -sf {input} {output}"


rule sumstats_link:
	input: lambda wildcards: config["sumstats"][wildcards.cohort]
	output: "sumstats/{cohort}.gz"
	shell: "ln -sf {input} {output}"

rule sumstats:
	input: "sumstats/ukb_mdd.md.eur.glm.logistic.gz",
		"sumstats/Depressio_FinnGen_R5_18032020.txt.gz"

rule daner_ukb:
	input: "sumstats/ukb_mdd.md.eur.glm.logistic.gz"
	output: "daner/daner_mdd_UKBB_MD.eur.glm.gz"
	shell: "Nca=$(zcat {input} | sed -n '2p' | awk '{{print $6/2}}');"
		"Nco=$(zcat {input} | sed -n '2p' | awk '{{print $6/2}}');"
		"echo \"CHR SNP BP A1 A2 FRQ_A_${{Nca}} FRQ_U_${{Nco}} INFO OR SE P\" > {output}"
              

rule daner:
	input: "daner/daner_MDD29.0515a_mds6.0316.gz",
		"daner/daner_GERA.euro.depress.0915a_mds5.id.gz",
		"daner/daner_mdd_decode_160211.gz",
		"daner/daner_mdd_genscot_1215a.gz",
		"daner/daner_mddGWAS_new_ipsych_170220.meta.gz",
		"daner/daner_mdd_23andMe_eur_v7.2.gz"
rule align:
	input: daner="daner/{cohort}.gz", ref="reference_info"
	output: "aligned/{cohort}.aligned.gz"
	script: "scripts/align.R"

rule refdir:
	output: "reference_info"
	shell: "impute_dirsub --refdir {config[refdir]} --reference_info --outname meta"

rule meta:
	input: "aligned/{cohort}.aligned.gz",
