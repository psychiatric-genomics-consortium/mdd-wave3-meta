# shared resources

## Download 1000 Genome Resources
rule resources_1kg_bed_download:
	input: HTTP.remote("https://www.dropbox.com/s/k9ptc4kep9hmvz5/1kg_phase1_all.tar.gz?dl=1", keep_local=False)
	output: "resources/1kg/1kg_phase1_all.tar.gz"
	shell: "cp {input} {output}"

rule resources_1kg_phase3_psam_download:
	input: HTTP.remote("https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1", keep_local=False)
	output: "resources/1kg/phase3_corrected.psam"
	shell: "cp {input} {output}"	
    
rule resources_1kg_phase3_pgen_download:
    input: HTTP.remote("https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1", keep_local=False)
    output: temp("resources/1kg/all_phase3.pgen.zst")
    shell: "cp {input} {output}"	
    
rule resources_1kg_phase3_pgen_decompress:
    input: "resources/1kg/all_phase3.pgen.zst"
    output: "resources/1kg/all_phase3.pgen"
    conda: "../envs/cojo.yaml"
    shell: "plink2 --zst-decompress {input} > {output}"

rule resources_1kg_phase3_pvar_download:
    input: HTTP.remote("https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1", keep_local=False)
    output: "resources/1kg/all_phase3.pvar.zst"
    shell: "cp {input} {output}"	
	
rule resources_1kg_bed:
	input: "resources/1kg/1kg_phase1_all.tar.gz"
	output: "resources/1kg/1kg_phase1_all.bed", "resources/1kg/1kg_phase1_all.bim", "resources/1kg/1kg_phase1_all.fam"
	shell: "tar --directory $(dirname {input}) -xzf {input}"
    
# find duplicate SNP IDs 
rule resource_1kg_dups:
    input: "resources/1kg/1kg_phase1_all.bim"
    output: "resources/1kg/1kg_phase1_all.dups"
    shell: "cat {input} | awk '{{print $2}}' | sort | uniq -d > {output}"
    
## Calculate ancestries frequencies for 1000G files
rule resources_1kg_phase3_frq:
    input: psam="resources/1kg/phase3_corrected.psam", pvar="resources/1kg/all_phase3.pvar.zst", pgen="resources/1kg/all_phase3.pgen"
    params: prefix="resources/1kg/phase3.{ancestries}", popname=lambda wildcards: wildcards.ancestries.upper()
    output: "resources/1kg/phase3.{ancestries}.afreq"
    conda: "../envs/cojo.yaml"
    shell: """plink2 --psam {input.psam} \
    --pvar {input.pvar} \
    --pgen {input.pgen} \
    --freq 'cols'='chrom,pos,ref,alt,reffreq,altfreq,nobs' \
    --keep-if SuperPop == {params.popname} \
    --rm-dup 'exclude-all' \
    --max-alleles 2 \
    --out {params.prefix}
    """
