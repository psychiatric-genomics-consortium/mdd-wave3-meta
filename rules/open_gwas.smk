# Query the IEU Open GWAS database

# Install the Open GWAS R library
rule open_gwas_install_ieugwasr:
    output: touch("resources/open_gwas/install_ieugwasr.done")
    conda: "../envs/open.yaml"
    shell: """
    Rscript -e "devtools::install_github('mrcieu/ieugwasr')"
    """
    
# Install the GWAS VCF R library
rule open_gwas_install_gwasvcf:
    output: touch("resources/open_gwas/install_gwasvcf.done")
    conda: "../envs/open.yaml"
    shell: """
    Rscript -e "devtools::install_github('mrcieu/gwasvcf')"
    """

# Do a lookup of all top SNPs 
# Before running this command, you'll need to fetch an authentication token using the
# command:
#  ieugwasr::get_access_token()
# in R
rule open_gwas_phewas_lookup:
    input: "docs/tables/meta_snps_full_eur.cojo.txt", ieugwas="resources/open_gwas/install_ieugwasr.done"
    output: phewas="results/open_gwas/phewas.txt", gwasinfo="results/open_gwas/gwasinfo.txt"
    conda: "../envs/open.yaml"
    script: "../scripts/open/phewas.R"
    
# fetch sumstats VCF files for GWAS matched by phewas query by dataset ID
rule open_gwas_fetch_dataset:
    input: dataset=HTTP.remote("https://gwas.mrcieu.ac.uk/files/{dataset}/{dataset}.vcf.gz"), phewas="results/open_gwas/phewas.txt", gwasinfo="results/open_gwas/gwasinfo.txt"
    output: temp("resources/open_gwas/datasets/{dataset}.vcf.gz")
    shell: "cp {input.dataset} {output}"
    
rule open_gwas_fetch_dataset_index:
    input: dataset=HTTP.remote("https://gwas.mrcieu.ac.uk/files/{dataset}/{dataset}.vcf.gz.tbi"), phewas="results/open_gwas/phewas.txt", gwasinfo="results/open_gwas/gwasinfo.txt"
    output: temp("resources/open_gwas/datasets/{dataset}.vcf.gz.tbi")
    shell: "cp {input.dataset} {output}"

# convert to text file
rule open_gwas_vcf2txt:
    input: vcf="resources/open_gwas/datasets/{dataset}.vcf.gz", tbi="resources/open_gwas/datasets/{dataset}.vcf.gz", hm3="resources/ldsc/w_hm3.snplist", gwasvcf_install="resources/open_gwas/install_gwasvcf.done"
    output: temp("resources/open_gwas/datasets/{dataset}.hm3.txt.gz")
    conda: "../envs/open.yaml"
    script: "../scripts/open/vcf2txt.R"

# munge for LDSC
rule open_gwas_munge:
    input: txt="resources/open_gwas/datasets/{dataset}.hm3.txt.gz", ldsc=rules.ldsc_install.output, hm3="resources/ldsc/w_hm3.snplist"
    output: sumstats="resources/open_gwas/datasets/{dataset}.hm3.sumstats.gz"
    params: prefix="resources/open_gwas/datasets/{dataset}.hm3"
    conda: "../envs/ldsc.yaml"
    shell: "resources/ldsc/ldsc/munge_sumstats.py --sumstats {input.txt} --out {params.prefix} --merge-alleles {input.hm3} --chunksize 500000"
    
# genetic correlation
rule open_gwas_ldsc_rg:
    input: mdd="results/distribution/daner_pgc_mdd_full_eur_hg19_v{version}.gz.ldsc.sumstats.gz",other="resources/open_gwas/datasets/{dataset}.hm3.sumstats.gz", w_ld=rules.ldsc_unzip_eur_w_ld_chr.output
    output: "results/open_gwas/ldsc/full_eur_v{version}_{dataset}.rg.log"
    params: prefix="results/open_gwas/ldsc/full_eur_v{version}_{dataset}.rg"
    conda: "../envs/ldsc.yaml"
    shell: "resources/ldsc/ldsc/ldsc.py --rg {input.mdd},{input.other} --ref-ld-chr {input.w_ld}/ --w-ld-chr {input.w_ld}/ --out {params.prefix}"

# return ids of GWAS datasets with SNP counts >= 2M
def open_gwas_parse_dataset_ids(gwasinfo_txt, min_snps=2e6):
    import pandas as pd
    import os.path
    if(os.path.exists(gwasinfo_txt)):
        gwasinfo = pd.read_table(gwasinfo_txt)
        return(gwasinfo[gwasinfo.nsnp >= min_snps].id.values)
    else:
        return([])
    
# list all IDs to fetch
rule open_gwas_rg_all_datasets:
    input: expand("results/open_gwas/ldsc/full_eur_v{version}_{dataset}.rg.log", version=analysis_version, dataset=open_gwas_parse_dataset_ids(rules.open_gwas_phewas_lookup.output.gwasinfo))

