# METAL meta-analysis with overlap correction

# Download and build METAL
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

# format sumstats for input into METAL
rule metal_sumstats_tbl:
	input: "results/sumstats/aligned/daner_{sumstats}.aligned.gz"
	output: "results/sumstats/tbl/{sumstats}.tbl"
	shell: """
	zcat {input} | awk -v OFS='\t' -v Neff=$(zcat {input} | head -n 1 | awk '{{print $6, $7}}' | sed "s/_/ /g" | awk '{{print (4*$3*$6)/($3 + $6)}}') '{{if(NR == 1) {{print $0, "FRQ_A", "Neff"}} else {{if(NF == 19) {{print $0, $6, 2*$19}} else {{print $0, $6, Neff}}}}}}' > {output}
	"""
	
rule metal_script_eur:
	input: tbls=expand("results/sumstats/tbl/mdd_{cohort}.eur.hg19.{release}.tbl", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur]), metal="resources/metal/METAL/build/bin/metal"
	params: process=lambda wildcards, input: '\n'.join(['PROCESS ' + tbl for tbl in input.tbls]), outfile="results/meta/metal/mdd_eur_v{version}_"
	output: "results/meta/metal/mdd_eur_v{version}.metal"
	run:
		metal = """
MARKER SNP
ALLELE A1 A2
PVALUE P
EFFECT log(OR)
WEIGHT Neff
CHROMOSOME CHR
POSITION BP
FREQ FRQ_A

SCHEME SAMPLESIZE
TRACKPOSITIONS ON
AVERAGEFREQ ON
OVERLAP ON

{process}

OUTFILE {outfile} .tbl
ANALYZE HETEROGENEITY
QUIT
		""".format(process=params.process, outfile=params.outfile)
		with open(output[0], 'w') as out:
			print(metal, file=out)
			
rule metal_eur:
	input: script="results/meta/metal/mdd_eur_v{version}.metal", metal="resources/metal/METAL/build/bin/metal"
	output: "results/meta/metal/mdd_eur_v{version}_1.tbl"
	shell: "{input.metal} {input.script}"
	
rule metal_eur_analyze:
	input: expand("results/meta/metal/mdd_eur_v{version}_1.tbl", version=analysis_version)
	
rule metal_daner:
	input: "results/meta/metal/mdd_eur_v{version}_1.tbl"
	output: "results/meta/metal/daner_mdd_eur_v{version}.metal.gz"
	shell: """
	cat {input} | awk -v OFS='\\t' '{{if(NR == 1) {{print "CHR", "SNP", "BP", "A1", "A2", "FRQ_A", "FRQ_U", "INFO", "BETA", "SE", "P"}} else {{print $1, $3, $2, toupper($4), toupper($5), $6, $6, 1, $9, 1, $11}}}}' | gzip -c > {output}
	"""