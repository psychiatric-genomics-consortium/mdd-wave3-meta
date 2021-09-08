
###################################
#                                 #
# Summary statistics distribution #
#                                 #
###################################


# Distribute results
# extensions and prefixes of Ricopili distribution output files
# daner_pgc_mdd_full_eur_hg19_v{version}.EXT
# Sumstats files
distribution_daner_ext = ["gz", "rp.gz", "gz.ldsc.sumstats.gz", "gz.p3.gz"]
# Clumped summary files
distribution_daner_p4_ext = ["gz.p4.clump.areator.sorted.1mhc", "gz.p4.clump.areator.sorted.1mhc.pvsorted", "gz.p4.clump.areator.sorted.1mhc.pvsorted.regs.txt", "gz.p4.clump.areator.sorted.1mhc.summary", "gz.p4.clump.areator.sorted.1mhc.xls", "het.gz.p4.clump.areator.sorted.1mhc", "het.gz.p4.clump.areator.sorted.1mhc.xls"]

# PREFIX.pgc_mdd_full_eur_hg19_v{version}.pdf
distribution_pdf_prefix = ["areas.fo", "areas", "manhattan.nog2", "manhattan.nog", "manhattan.v2", "qq"]

# PREFIX.pgc_mdd_full_eur_hg19_v{version}.het.pdf
distribution_het_pdf_prefix = ["manhattan.v2", "qq"]

# basic.pgc_mdd_full_eur_hg19_v{version}.EXT
distribution_basic_ext = ["num.xls"]

# Distribute meta analysis files
rule distribute_meta:
	input: "results/meta/distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{file}"
	output: DBox_dist.remote("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{file}")
	shell: "cp {input} {output}"
	
	
# share files locally
distribute_local_path=config["remote"]["distribution"]["lisa"] if "lisa" in config["remote"]["distribution"] else os.path.expanduser('~')
rule distribute_local:
	input: "results/meta/distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{file}"
	output: expand("{local_path}/mdd3/distribution/pgc_mdd_{{cohorts}}_eur_hg19_v{{version}}/{{file}}", local_path=distribute_local_path)
	shell: "cp {input} {output}"

##
## Distribution for PGC analysts
##

# list all files to be uploaded to Dropbox
rule DBox_dist_analyst:
	input: DBox_dist.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/daner_pgc_mdd_{cohorts}_eur_hg19_v{version}.{ext}", version=analysis_version, cohorts=cohorts_analyst, ext=distribution_daner_ext))
	DBox_dist.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/daner_pgc_mdd_{cohorts}_eur_hg19_v{version}.{ext}", version=analysis_version, cohorts=cohorts_analyst, ext=distribution_daner_p4_ext)),
	DBox_dist.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{prefix}.pgc_mdd_{cohorts}_eur_hg19_v{version}.pdf", version=analysis_version, cohorts=cohorts_analyst, prefix=distribution_pdf_prefix)),
	DBox_dist.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{prefix}.pgc_mdd_{cohorts}_eur_hg19_v{version}.het.pdf", version=analysis_version, cohorts=cohorts_analyst, prefix=distribution_het_pdf_prefix)),
	DBox_dist.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/basic.pgc_mdd_{cohorts}_eur_hg19_v{version}.num.xls", version=analysis_version, cohorts=cohorts_analyst, ext=distribution_basic_ext))
	
# list all files to be shared locally (e.g., on LISA)
rule local_dist_analyst:
	input: expand("{local_path}/mdd3/distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/daner_pgc_mdd_{cohorts}_eur_hg19_v{version}.{ext}", local_path=distribute_local_path, version=analysis_version, cohorts=cohorts_analyst, ext=distribution_daner_ext),
	expand("{local_path}/mdd3/distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/daner_pgc_mdd_{cohorts}_eur_hg19_v{version}.{ext}", local_path=distribute_local_path, version=analysis_version, cohorts=cohorts_analyst, ext=distribution_daner_p4_ext),
	expand("{local_path}/mdd3/distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{prefix}.pgc_mdd_{cohorts}_eur_hg19_v{version}.pdf", local_path=distribute_local_path, version=analysis_version, cohorts=cohorts_analyst, prefix=distribution_pdf_prefix),
	expand("{local_path}/mdd3/distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{prefix}.pgc_mdd_{cohorts}_eur_hg19_v{version}.het.pdf", local_path=distribute_local_path, version=analysis_version, cohorts=cohorts_analyst, prefix=distribution_het_pdf_prefix),
	expand("{local_path}/mdd3/distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/basic.pgc_mdd_{cohorts}_eur_hg19_v{version}.num.xls", local_path=distribute_local_path, version=analysis_version, cohorts=cohorts_analyst, ext=distribution_basic_ext)

# Download daner sumstats for downstream analysis
# Look at config file to determine whether to fetch locally on LISA or remotely
# from Dropbox share
rule redistribute_daner:
	input: lambda wildcards: expand("{local_path}/mdd3/distribution/{analysis}/daner_{analysis}.{{ext}}.gz", local_path=distribute_local_path, analysis=wildcards.analysis) if "lisa" in config["remote"]["distribution"] else
		DBox_dist.remote(expand("distribution/{analysis}/daner_{analysis}.{{ext}}.gz", analysis=wildcards.analysis))
	output: "results/distribution/daner_{analysis}.{ext}.gz"
	shell: "cp {input} {output}"

rule downstream_full:
	input: expand("results/distribution/daner_pgc_mdd_full_eur_hg19_v{version}.rp.gz", version=analysis_version)

rule downstream_noUKBB:
	input: expand("results/distribution/daner_pgc_mdd_noUKBB_eur_hg19_v{version}.rp.gz", version=analysis_version)
	
rule downstream_noALSPAC:
	input: expand("results/distribution/daner_pgc_mdd_noALSPAC_eur_hg19_v{version}.rp.gz", version=analysis_version)

rule downstream_no23andMe:
	input: expand("results/distribution/daner_pgc_mdd_no23andMe_eur_hg19_v{version}.rp.gz", version=analysis_version)

# Download tables and figures
rule redistribute_figtabs:
	input: lambda wildcards: expand("{local_path}/mdd3/distribution/{analysis}_v{version}/{prefix}.{analysis}_v{version}.{ext}", local_path=distribute_local_path, analysis=wildcards.analysis, version=wildcards.version, ext=wildcards.ext) if "lisa" in config["remote"]["distribution"] else DBox_dist.remote("distribution/{analysis}_v{version}/{prefix}.{analysis}_v{version}.{ext}")
	output: "results/distribution/{prefix}.{analysis}_v{version}.{ext}"
	shell: "cp {input} {output}"

rule redistribute_danerxls:
	input: lambda wildcards: expand("{local_path}/mdd3/distribution/{analysis}_v{version}/daner_{analysis}_v{version}.xls", local_path=distribute_local_path, analysis=wildcards.analysis, version=wildcards.version, ext=wildcards.ext) if "lisa" in config["remote"]["distribution"] else DBox_dist.remote("distribution/{analysis}_v{version}/daner_{analysis}_v{version}.xls")
	output: "results/distribution/daner_{analysis}_v{version}.xls"
	shell: "cp {input} {output}"
	
rule redistribute_danerext:
	input: lambda wildcards: expand("{local_path}/mdd3/distribution/{analysis}_v{version}/daner_{analysis}_v{version}.gz.{ext}", local_path=distribute_local_path, analysis=wildcards.analysis, version=wildcards.version, ext=wildcards.ext) if "lisa" in config["remote"]["distribution"] else DBox_dist.remote("distribution/{analysis}_v{version}/daner_{analysis}_v{version}.gz.{ext}")
	output: "results/distribution/daner_{analysis}_v{version}.gz.{ext}"
	shell: "cp {input} {output}"
	
# alternative versions (like METACARPA)
rule redistribute_daner_alt_ext:
	input: lambda wildcards: expand("{local_path}/mdd3/distribution/{analysis}_v{version}/daner_{analysis}_v{version}.{{alt}}.gz.{ext}", local_path=distribute_local_path, analysis=wildcards.analysis, version=wildcards.version, ext=wildcards.ext) if "lisa" in config["remote"]["distribution"] else DBox_dist.remote("distribution/{analysis}_v{version}/daner_{analysis}_v{version}.{alt}.gz.{ext}")
	output: "results/distribution/daner_{analysis}_v{version}.{alt}.gz.{ext}"
	shell: "cp {input} {output}"

# download most recent manhattan plot
rule manhattan_full:
	input: expand("results/distribution/manhattan.nog2.pgc_mdd_full_eur_hg19_v{version}.pdf", version=analysis_version)


##
## Distribution for public
##

ruleorder: distribute_public > distribute_meta

rule distribute_public:
	input: "results/meta/distribution/pgc_mdd_no23andMe_eur_hg19_v{version}/{file}"
	output: DBox_dist_public.remote("distribution/pgc_mdd_no23andMe_eur_hg19_v{version}/{file}")
	shell: "cp {input} {output}"

# list all files to be uploaded to Dropbox
rule DBox_dist_public:
	input: DBox_dist_public.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/daner_pgc_mdd_{cohorts}_eur_hg19_v{version}.{ext}", version=analysis_version, cohorts=cohorts_public, ext=distribution_daner_ext)), DBox_dist_public.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{prefix}.pgc_mdd_{cohorts}_eur_hg19_v{version}.pdf", version=analysis_version, cohorts=cohorts_public, prefix=distribution_pdf_prefix)), DBox_dist_public.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/{prefix}.pgc_mdd_{cohorts}_eur_hg19_v{version}.het.pdf", version=analysis_version, cohorts=cohorts_public, prefix=distribution_het_pdf_prefix)), DBox_dist_public.remote(expand("distribution/pgc_mdd_{cohorts}_eur_hg19_v{version}/basic.pgc_mdd_{cohorts}_eur_hg19_v{version}.num.xls", version=analysis_version, cohorts=cohorts_public, ext=distribution_basic_ext))


