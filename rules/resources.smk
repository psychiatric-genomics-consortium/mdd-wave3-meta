# shared resources

rule resources_1kg_bed_download:
	input: HTTP.remote("https://www.dropbox.com/s/k9ptc4kep9hmvz5/1kg_phase1_all.tar.gz?dl=1", keep_local=True)
	output: "resources/1kg/1kg_phase1_all.tar.gz"
	shell: "cp {input} {output}"

rule resources_1kg_psam_download:
	input: HTTP.remote("https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1", keep_local=True)
	output: "resources/1kg/phase3_corrected.psam"
	shell: "cp {input} {output}"	
	
rule resources_1kg_bed:
	input: "resources/1kg/1kg_phase1_all.tar.gz"
	output: "resources/1kg/1kg_phase1_all.bed", "resources/1kg/1kg_phase1_all.bim", "resources/1kg/1kg_phase1_all.fam"
	shell: "tar --directory $(dirname {input}) -xzf  {input}"