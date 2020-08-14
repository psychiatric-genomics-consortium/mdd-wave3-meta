# Store GWAS data in VCF format https://github.com/MRCIEU/gwas-vcf-specification

# install the vcfgwas library
rule install_github_vcfgwas:
    output: "logs/vcf/install_github.log"
    conda: "../envs/vcf.yaml"
    shell: "Rscript -e 'devtools::install_github(\"mrcieu/gwasvcf\")' > {output}"