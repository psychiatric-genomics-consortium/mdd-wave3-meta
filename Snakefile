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
		"sumstats/Depressio_FinnGen_R5_18032020.txt.gz",
		"sumstats/ALSPAC_FILTERED_ALL_12082019.gz"

rule daner_ukb:
	input: "sumstats/ukb_mdd.md.eur.glm.logistic.gz"
	output: "daner/daner_mdd_UKBB_MD.eur.glm.gz"
	shell: "Nca=$(zcat {input} | sed -n '2p' | awk '{{print $6/2}}');"
		"Nco=$(zcat {input} | sed -n '2p' | awk '{{print $7/2}}');"
   		"echo \"CHR SNP BP A1 A2 FRQ_A_${{Nca}} FRQ_U_${{Nco}} INFO OR SE P\" > daner/daner_mdd_UKBB_MD.eur.glm;"
   		"zcat {input} | tail -n +2 | awk '{{print $1, $3, $2, $4, $5, $8, $9, $10, $13, $14, $15}}' >> daner/daner_mdd_UKBB_MD.eur.glm;"
		"gzip --verbose daner/daner_mdd_UKBB_MD.eur.glm"
              
rule daner_finn:
	input: "sumstats/Depressio_FinnGen_R5_18032020.txt.gz"
        output: "liftover/daner_mdd_FinnGen_R5.gz"
	shell: "echo \"CHR SNP BP A1 A2 FRQ_A_23424 FRQ_U_192220 INFO OR SE P\" > daner/daner_mdd_FinnGen_R5;"
		"zcat {input} | tail -n +2 | awk '{{print $1, $5, $2, $4, $3, $11, $12, $18, exp($8), $9, $7}}' >> daner/daner_mdd_FinnGen_R5;"
		"gzip --verbose daner/daner_mdd_FinnGen_R5"

rule liftover_finn:
	input: "liftover/daner_mdd_FinnGen_R5.gz"
	output: "daner/daner_mdd_FinnGen_R5.gz"
	script: "scripts/liftover.R"

rule daner:
	input: "daner/daner_MDD29.0120a.rmUKBB.gz",
   		"daner/daner_GERA.euro.depress.0915a_mds5.id.gz",
		"daner/daner_mdd_decode_160211.gz",
   		"daner/daner_mdd_genscot_1215a.gz",
   		"daner/daner_mddGWAS_new_ipsych_170220.meta.gz",
   		"daner/daner_mdd_23andMe_eur_v7.2.gz"

rule align:
	input: daner="daner/{cohort}.gz", ref="meta/reference_info"
	output: "aligned/{cohort}.aligned.gz"
	script: "scripts/align.R"

rule refdir:
	output: "meta/reference_info"
	shell: "cd meta; impute_dirsub --refdir {config[refdir]} --reference_info --outname meta"

rule meta:
	input: "aligned/{cohort}.gz"
	output: "meta/{cohort}.gz"
	shell: "ln -sf $(readlink -f {input}) {output}"

rule dataset_eur:
	input: "meta/daner_MDD29.0120a.rmUKBB.aligned.gz", "meta/daner_GERA.euro.depress.0915a_mds5.id.aligned.gz", "meta/daner_mdd_decode_160211.aligned.gz", "meta/daner_mdd_genscot_1215a.aligned.gz", "meta/daner_mddGWAS_new_ipsych_170220.meta.aligned.gz", "meta/daner_mdd_23andMe_eur_v7.2.aligned.gz"
        output: "meta/dataset_eur"
	shell: "for daner in {input}; do echo $(basename $daner) >> {output}; done"

rule postimp_eur:
	input: dataset="meta/dataset_eur", ref="meta/reference_info"
	output: "meta/j._pi_pgc_mdd_meta_eur_v3.00_2020.07.id"
	shell: "cd meta; postimp_navi --result {input.dataset} --popname eur --out pgc_mdd3_meta_eur_v3.00_2020.07"
