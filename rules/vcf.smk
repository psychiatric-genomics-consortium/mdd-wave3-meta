# Store GWAS data in VCF format https://github.com/MRCIEU/gwas-vcf-specification (vcf.gz)
# and PGC VCF-like format (pgc.gz)

# install gwas2vcf
rule vcf_install_gwas2vcf:
    output: directory("resources/vcf/vendor/gwas2vcf")
    conda: "../envs/vcf.yaml"
    shell: "pip install git+git://github.com/bioinformed/vgraph@v1.4.0#egg=vgraph; git clone git@github.com:MRCIEU/gwas2vcf.git {output}"

# fasta reference files
rule vcf_fasta_grch37:
    input: HTTP.remote("http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta")
	output: "resources/fasta/human_grch37.fasta"
	shell: "cp {input} {output}"

rule vcf_fasta_grch37_fai:
    input: HTTP.remote("http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.fai")
	output: "resources/fasta/human_grch37.fasta.fai"
	shell: "cp {input} {output}"

rule vcf_fasta_grch38:
    input: HTTP.remote("https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta")
	output: "resources/fasta/human_grch38.fasta"
	shell: "cp {input} {output}"

rule vcf_fasta_grch38_fai:
    input: HTTP.remote("https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai")
	output: "resources/fasta/human_grch38.fasta.fai"
	shell: "cp {input} {output}"


# dbSNP reference files
# GRCh37/hg19/b37
rule vcf_dbsnp_grch37:
    input: HTTP.remote("http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.b37.vcf.gz")
    output: "resources/dbsnp/human_grch37.dbsnp.v153.vcf.gz"
    shell: "cp {input} {output}"

rule vcf_dbsnp_grch37_tbi:
    input: HTTP.remote("http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.b37.vcf.gz.tbi")
    output: "resources/dbsnp/human_grch37.dbsnp.v153.vcf.gz.tbi"
    shell: "cp {input} {output}"

rule vcf_dbsnp_grch38:
    input: HTTP.remote("http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.gz")
    output: "resources/dbsnp/human_grch38.dbsnp.vcf.gz"
    shell: "cp {input} {output}"

rule vcf_dbsnp_grch38_tbi:
    input: HTTP.remote("http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.gz.tbi")
    output: "resources/dbsnp/human_grch38.dbsnp.v153.vcf.gz.tbi"
    shell: "cp {input} {output}"

# Convert OR to Log-Odds
# convert CHR "23" to "X" to match fasta files
rule vcf_logOR:
    input: "results/distribution/daner_{sumstats}.rp.gz"
    output: "results/vcf/beta/{sumstats}.gz"
    shell: "gunzip -c {input} | tail -n +2 | awk '{{if($1 == 23) ($1 = \"X\"); print $1, $3, $4, $5, log($9), $10, $11, $2, $7, $8, $17, $18}}' | gzip -c > {output}"

# Parameter file for VCF conversion (json)
builds = {"19": "GRCh37", "38": "GRCh38"}
rule vcf_daner2vcf_json:
    input: daner="results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.rp.gz"
    output: vcf="results/vcf/gwas/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.json"
    log: "logs/vcf/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.json.log"
    run:
        with gzip.open(input.daner, 'r') as daner:
            headers = daner.readline().split()
            cohort_cases = headers[5].decode().split('_')[2]
            cohort_controls = headers[6].decode().split('_')[2]
        with open(output.vcf, 'w') as out:
            json.dump({"chr_col": 0,
                        "pos_col": 1,
                        "ea_col": 2,
                        "oa_col": 3,
                        "beta_col": 4,
                        "se_col": 5,
                        "pval_col": 6,
                        "snp_col": 7,
                        "eaf_col": 8,
                        "imp_info_col": 9,
                        "ncase_col": 10,
                        "ncontrol_col": 11,
                        "delimiter": " ",
                        "header": False,
                        "build": "{build}".format(build=builds[wildcards.hg]),
                        "cohort_cases": cohort_cases,
                        "cohort_controls": cohort_controls,
                        "id": "pgc-mdd2024-{cohorts}-{pops}".format(cohorts=wildcards.cohorts, pops=wildcards.ancestries)},
                     out)


# Convert daner to vcf
rule vcf_daner2vcf:
    input: beta="results/vcf/beta/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.gz", json="results/vcf/gwas/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.json", fasta=lambda wildcards: expand("resources/fasta/human_{build}.fasta", build=builds[wildcards.hg].lower()), fai=lambda wildcards: expand("resources/fasta/human_{build}.fasta.fai", build=builds[wildcards.hg].lower()), dbsnp=lambda wildcards: expand("resources/dbsnp/human_{build}.dbsnp.v153.vcf.gz", build=builds[wildcards.hg].lower()), dbsnp_tbi=lambda wildcards: expand("resources/dbsnp/human_{build}.dbsnp.v153.vcf.gz.tbi", build=builds[wildcards.hg].lower()), gwas2vcf=rules.vcf_install_gwas2vcf.output
    output: "results/vcf/gwas/pgc-mdd{year}_{cohorts}_{ancestries}_hg{hg}_v{analysis}.vcf.gz"
    log: "logs/vcf/pgc-mdd{year}_{cohorts}_{ancestries}_hg{hg}_v{analysis}.log"
    conda: "../envs/vcf.yaml"
    shell: "python resources/vcf/vendor/gwas2vcf/main.py --out {output} --data {input.beta} --json {input.json} --ref {input.fasta} --dbsnp {input.dbsnp} > {log}"
	
# Annotation for effective sample size
rule vcf_neff_annotation:
	input: daner="results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.rp.gz"
	output: annot="results/vcf/annot/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.ne.gz", annot_tbi="results/vcf/annot/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.ne.gz.tbi"
	shell: """
	gunzip -c {input.daner} | awk 'OFS="\\t" {{if(NR == 1) {{print "#CHROM", "POS", "NE"}} else {{if($1 == "23") $1 == "X"; print $1, $3, 2*$19}}}}' | bgzip -c > {output.annot}
	tabix -s1 -b2 -e2 {output.annot}
	"""
	
# Add NE annotation to VCF FORMAT for the GWAS sample
rule vcf_neff_annotate:
	input: annot="results/vcf/annot/pgc_mdd_{cohorts}_{ancestries}_hg19_v{analysis}.ne.gz", annot_tbi="results/vcf/annot/pgc_mdd_{cohorts}_{ancestries}_hg19_v{analysis}.ne.gz.tbi", vcf="results/vcf/gwas/pgc-mdd{year}_{cohorts}_{ancestries}_hg19_v{analysis}.vcf.gz"
	output: vcf="results/vcf/pgc-mdd{year}_{cohorts}_{ancestries}_v{analysis}.vcf.gz"
	params: id=lambda wildcards: "pgc-mdd{year}-{cohorts}-{pops}".format(cohorts=wildcards.cohorts, pops=wildcards.ancestries, year=wildcards.year)
	shell: """
	echo -e '##FORMAT=<ID=NE,Number=1,Type=Integer,Description="Effective sample size used to estimate genetic effect">' >> {input.annot}.hdr.txt
	bcftools annotate -s {params.id} -a {input.annot} -h {input.annot}.hdr.txt -c CHROM,POS,FORMAT/NE -o  {output.vcf} {input.vcf}
	"""

# List sumstats to convert to VCF
rule vcf:
    input: expand("results/vcf/pgc-mdd{year}_{cohorts}_{ancestries}_v{analysis}.vcf.gz", year=2023, cohorts=['full', 'no23andMe', 'Clin', 'EHR', 'Quest', 'SelfRep'], ancestries=['eur'], analysis=analysis_version)
    
# VCF-like PGC sumstats file
# Pull in daner file for sumstats, fai file for genome build info, basic.num file
# for genotyped and sumstats cohorts with case/control/allele counts
rule vcf_daner2pgc:
    input: daner="results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{wave}.{geno}.{sums}.{rev}.neff.gz",fasta_fai="resources/fasta/human_grch37.fasta.fai", cohorts="docs/tables/cohorts/{cohorts}_{ancestries}_v{wave}.{geno}.{sums}.{rev}.txt",  cff="CITATION.cff", header_template="scripts/vcf/pgc.glue"
    conda: "../envs/vcf.yaml"
    output: temp("results/vcf/pgc-mdd{year}_{cohorts}_{ancestries}_v{wave}-{geno}-{sums}-{rev}.tsv")
    script: "../scripts/vcf/pgc.R"

rule vcf_pgc_gz:
    input: "results/vcf/pgc-mdd{sumstats}.tsv"
    output: "results/vcf/pgc-mdd{sumstats}.tsv.gz"
    shell: "gzip -c {input} > {output}"  
    
rule vcf_pgc:
    input: expand("results/vcf/pgc-mdd{year}_{analysis}.tsv.gz", year=2024, analysis=['full_div_v3-49-46-01', 'full_eur_v3-49-24-11', 'no23andMe_div_v3-49-46-01', 'no23andMe_eur_v3-49-24-11', 'no23andMe-noUKBB_eur_v3-49-24-11', 'Clin_eur_v3-49-24-11', 'EHR_eur_v3-49-24-11', 'Quest_eur_v3-49-24-11', 'SelfRep_eur_v3-49-24-11', 'top10k_eur_v3-49-24-11'])