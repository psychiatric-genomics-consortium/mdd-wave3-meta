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