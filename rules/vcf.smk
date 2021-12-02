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


# Conversion of final sumstats
vcf_sumstats_gz, = glob_wildcards("results/distribution/daner_{sumstats}.gz")

# Convert OR to Log-Odds
rule vcf_logOR:
    input: "results/distribution/daner_{sumstats}.gz"
    output: "results/vcf/beta/{sumstats}.gz"
    shell: "gunzip -c {input} | tail -n +2 | awk '{{print $1, $3, $4, $5, log($9), $10, $11, $2, $7, $8, $17, $18}}' | gzip -c > {output}"

# Parameter file for VCF conversion (json)
builds = {"19": "GRCh37", "38": "GRCh38"}
rule vcf_daner2vcf_json:
    input: daner="results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.gz"
    output: vcf="results/vcf/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.json"
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
                        "id": "pgc_mdd_{cohorts}_{pops}_v{analysis}".format(analysis=wildcards.analysis, cohorts=wildcards.cohorts, pops=wildcards.ancestries.upper())},
                     out)



rule vcf_daner2vcf:
    input: beta="results/vcf/beta/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.gz", json="results/vcf/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.json", fasta=lambda wildcards: expand("resources/fasta/human_{build}.fasta", build=builds[wildcards.hg].lower()), fai=lambda wildcards: expand("resources/fasta/human_{build}.fasta.fai", build=builds[wildcards.hg].lower()), gwas2vcf=rules.vcf_install_gwas2vcf.output
    output: "results/vcf/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.vcf.gz"
    log: "logs/vcf/pgc_mdd_{cohorts}_{ancestries}_hg{hg}_v{analysis}.log"
    conda: "../envs/vcf.yaml"
    shell: "python resources/vendor/gwas2vcf/main.py --out {output} --data {input.beta} --json {input.json} --ref {input.fasta} > {log}"

rule vcf:
    input: expand("results/vcf/{sumstats}.vcf.gz", sumstats=vcf_sumstats_gz)
    
# VCF-like PGC sumstats file
rule vcf_daner2pgc:
    input: daner="results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{analysis}.neff.gz", fasta_fai="resources/fasta/human_grch37.fasta.fai", genotype_cohorts="docs/tables/cohorts_mdd.eur.txt", sumstats_cohorts="docs/tables/cohorts.eur.txt", cff="CITATION.cff", header_template="scripts/vcf/pgc.glue"
    conda: "../envs/vcf.yaml"
    output: "results/vcf/pgc-mdd{year}-{cohorts}-{ancestries}-v{analysis}.pgc"
    script: "../scripts/vcf/pgc.R"

rule vcf_pgc_gz:
    input: "results/vcf/pgc-mdd{sumstats}.pgc"
    output: "results/vcf/pgc-mdd{sumstats}.pgc.gz"
    shell: "gzip -c {input} > {output}"  
    
rule vcf_pgc:
    input: expand("results/vcf/pgc-mdd{year}-{cohorts}-{ancestries}-v{analysis}.pgc.gz", year=2022, cohorts=['full', 'no23andMe'], ancestries=['eur'], analysis=analysis_version)