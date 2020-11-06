# install LDSC
rule ldsc_install_hm3_bz:
	input: HTTP.remote("https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2")
	output: "resources/ldsc/w_hm3.snplist.bz2"
	shell: "cp {input} {output}"
	
rule meta_install_hm3:
	input: "resources/ldsc/w_hm3.snplist.bz2"
	output: "resources/ldsc/w_hm3.snplist"
	shell: "bunzip2 {input}"

	