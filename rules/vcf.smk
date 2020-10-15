# Store GWAS data in VCF format https://github.com/MRCIEU/gwas-vcf-specification

# install the vcfgwas library
rule vcf_install_github_vcfgwas:
    output: "resources/vcf/vendor/r-gwasvcf"
    log: "logs/vcf/install_github.log"
    conda: "../envs/vcf.yaml"
    shell: "Rscript -e 'devtools::install_github(\"mrcieu/gwasvcf\", upgrade=\"never\")' 2>&1 > {output}"

# install gwas2vcf
rule vcf_install_gwas2vcf:
    output: directory("resources/vcf/vendor/gwas2vcf")
    conda: "../envs/vcf.yaml"
    shell: "pip install git+git://github.com/bioinformed/vgraph@v1.4.0#egg=vgraph; git clone git@github.com:MRCIEU/gwas2vcf.git {output}"

# fasta reference files
rule vcf_fasta_grch37:
    input: HTTP.remote("http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta")
	output: "resources/fasta/human_g1k_v37.fasta"
	shell: "cp {input} {output}"

rule vcf_fasta_grch37_fai:
    input: HTTP.remote("http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.fai")
	output: "resources/fasta/human_g1k_v37.fasta.fai"
	shell: "cp {input} {output}"

rule vcf_fasta_grch38:
    input: HTTP.remote("https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta")
	output: "resources/fasta/Homo_sapiens_assembly38.fasta"
	shell: "cp {input} {output}"

rule vcf_fasta_grch38_fai:
    input: HTTP.remote("https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai")
	output: "resources/fasta/Homo_sapiens_assembly38.fasta.fai"
	shell: "cp {input} {output}"


vcf_sumstats_gz, = glob_wildcards("results/distribution/daner_{sumstats}.gz")

# Convert OR to Log-Odds
rule vcf_logOR:
    input: "results/distribution/daner_{sumstats}.gz"
    output: "results/vcf/beta/{sumstats}.gz"
    shell: "gunzip -c {input} | tail -n +2 | awk '{{print $1, $3, $4, $5, log($9), $10, $11, $2, $7, $8, $17, $18}}' | gzip -c > {output}"

# Parameter file for VCF conversion (json)
rule vcf_daner2vcf_json:
    input: daner="results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz"
    output: vcf="results/vcf/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.json"
    log: "logs/vcf/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.json.log"
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
                        "build": "GRCh37",
                        "cohort_cases": cohort_cases,
                        "cohort_controls": cohort_controls,
                        "id": "pgc_mdd_{cohorts}_{pops}_v{version}".format(version=wildcards.version, cohorts=wildcards.cohorts, pops=wildcards.ancestries.upper())},
                     out)



rule vcf_daner2vcf:
    input: beta="results/vcf/beta/{sumstats}.gz", json="results/vcf/{sumstats}.json", fasta="resources/fasta/human_g1k_v37.fasta", fai="resources/fasta/human_g1k_v37.fasta.fai", gwas2vcf=rules.vcf_install_gwas2vcf.output
    output: "results/vcf/{sumstats}.vcf.gz"
    log: "logs/vcf/{sumstats}.log"
    conda: "../envs/vcf.yaml"
    shell: "python resources/vendor/gwas2vcf/main.py --out {output} --data {input.beta} --json {input.json} --ref {input.fasta} > {log}"

rule vcf:
    input: expand("results/vcf/{sumstats}.vcf.gz", sumstats=vcf_sumstats_gz)