# install LDSC
rule ldsc_fetch_hm3_bz:
	input: HTTP.remote("https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2")
	output: "resources/ldsc/w_hm3.snplist.bz2"
	shell: "cp {input} {output}"
	
rule ldsc_unzip_hm3:
	input: "resources/ldsc/w_hm3.snplist.bz2"
	output: "resources/ldsc/w_hm3.snplist"
	shell: "bunzip2 {input}"

rule ldsc_fetch_eur_w_ld_chr_bz:
	input: HTTP.remote("https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2")
	output: "resources/ldsc/eur_w_ld_chr.tar.bz2"
	shell: "cp {input} {output}"

rule ldsc_unzip_eur_w_ld_chr:
	input: "resources/ldsc/eur_w_ld_chr.tar.bz2"
	output: directory("resources/ldsc/eur_w_ld_chr")
	shell: "tar -jxvf {input} -C $(dirname {output})"

rule ldsc_install:
	input: "resources/ldsc/w_hm3.snplist", rules.ldsc_unzip_eur_w_ld_chr.output
	output: directory("resources/ldsc/ldsc")
	shell: "git clone https://github.com/bulik/ldsc.git {output}"
	
rule ldsc_munge:
	input: sumstats="results/distribution/{cohort}.gz", hm3="resources/ldsc/w_hm3.snplist", ldsc=rules.ldsc_install.output
	params:
		prefix="results/ldsc/h2/{cohort}"
	conda: "../envs/ldsc.yaml"
	output: "results/ldsc/munged/{cohort}.sumstats.gz"
	shell: "resources/ldsc/ldsc/munge_sumstats.py --sumstats {input.sumstats} --daner --out {params.prefix} --merge-alleles {input.hm3} --chunksize 500000"	
	
rule ldsc_h2:
	input: sumstats="results/ldsc/munged/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.sumstats.gz", w_ld=rules.ldsc_unzip_eur_w_ld_chr.output
	params:
		prefix="results/ldsc/h2/pgc_mdd_{cohorts}_{ancestries}_v{version}"
	conda: "../envs/ldsc.yaml"
	output: "results/ldsc/h2/pgc_mdd_{cohorts}_{ancestries}_v{version}.log"
	shell: "resources/ldsc/ldsc/ldsc.py --h2 {input.sumstats} --ref-ld-chr {input.w_ld}/ --w-ld-chr {input.w_ld}/ --out {params.prefix}"

