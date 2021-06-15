# METAL meta-analysis with overlap correction

rule metal_install:
	input: HTTP.remote("http://csg.sph.umich.edu/abecasis/Metal/download/Linux-metal.tar.gz")
	output: dir=directory("resources/metal"), metal="resources/metal/generic-metal/metal"
	shell: "tar -C {output.dir} -xzf {input}"
	
rule metal_github_download:
	output: directory("resources/metal/METAL")
	shell: "git clone https://github.com/statgen/METAL {output}"
	
rule metal_build:
	output: "resources/metal/METAL/build/bin/metal"
	conda: "../envs/meta.yaml"
	shell: """
	cd resources/metal/METAL
	mkdir -p build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..
	make
	make install
	"""
	
	
rule metal_eur:
	input: sumstats=expand("results/sumstats/aligned/daner_mdd_{cohort}.eur.hg19.{release}.aligned.gz", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur]), metal="resources/metal/METAL/build/bin/metal"
	output: "results/meta/metal/mdd_{version}.table"
	shell: """
	{input.metal}
	
	MARKER	SNP
	ALLELE	A1 A2
	PVALUE	P
	EFFECT	log(OR)
	
	SCHEME SAMPLESIZE
	
	
	"""
	