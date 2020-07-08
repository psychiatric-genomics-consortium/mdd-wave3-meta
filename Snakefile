configfile: "config.yaml"

rule link:
	input: lambda wildcards: config["daner"][wildcards.cohort]
        output: "daner/{cohort}.gz"
	shell: "ln -sf {input} {output}"

rule daner:
	input: "daner/daner_MDD29.0515a_mds6.0316.gz",
               "daner/daner_GERA.euro.depress.0915a_mds5.id.gz",
	       "daner/daner_mdd_decode_160211.gz",
               "daner/daner_mdd_genscot_1215a.gz",
               "daner/daner_mddGWAS_new_ipsych_170220.meta.gz",
               "daner/daner_mdd_23andMe_eur_v7.2.gz"

