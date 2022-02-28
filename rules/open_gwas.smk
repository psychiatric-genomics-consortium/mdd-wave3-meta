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
# rule open_gwas_fetch_dataset:
#     input: dataset=HTTP.remote("https://gwas.mrcieu.ac.uk/files/{dataset}/{dataset}.vcf.gz"), phewas="results/open_gwas/phewas.txt", gwasinfo="results/open_gwas/gwasinfo.txt"
#     output: temp("resources/open_gwas/datasets/{dataset}.vcf.gz")
#     shell: "cp {input.dataset} {output}"
# rule open_gwas_fetch_dataset_index:
#     input: dataset=HTTP.remote("https://gwas.mrcieu.ac.uk/files/{dataset}/{dataset}.vcf.gz.tbi"), phewas="results/open_gwas/phewas.txt", gwasinfo="results/open_gwas/gwasinfo.txt"
#     output: temp("resources/open_gwas/datasets/{dataset}.vcf.gz.tbi")
#     shell: "cp {input.dataset} {output}"    
    
# use curl directly rather than the HTTP remote downloader so that tasks can be
# run in parallel    
rule open_gwas_fetch_dataset:
    output: temp("resources/open_gwas/datasets/{dataset}.vcf.gz")
    shell: "curl https://gwas.mrcieu.ac.uk/files/{wildcards.dataset}/{wildcards.dataset}.vcf.gz > {output}"
  
rule open_gwas_fetch_dataset_index:
    output: temp("resources/open_gwas/datasets/{dataset}.vcf.gz.tbi")
    shell: "curl https://gwas.mrcieu.ac.uk/files/{wildcards.dataset}/{wildcards.dataset}.vcf.gz.tbi > {output}"

# convert to text file
rule open_gwas_vcf2txt:
    input: vcf="resources/open_gwas/datasets/{dataset}.vcf.gz", tbi="resources/open_gwas/datasets/{dataset}.vcf.gz.tbi", hm3="resources/ldsc/w_hm3.snplist", gwasvcf_install="resources/open_gwas/install_gwasvcf.done", gwasinfo="results/open_gwas/gwasinfo.txt"
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
    input: mdd="results/ldsc/munged/pgc_mdd_{cohorts}_eur_hg19_v{version}.sumstats.gz",other="resources/open_gwas/datasets/{dataset}.hm3.sumstats.gz", w_ld=rules.ldsc_unzip_eur_w_ld_chr.output
    output: "results/open_gwas/ldsc/{cohorts}_eur_v{version}_{dataset}.rg.log"
    params: prefix="results/open_gwas/ldsc/{cohorts}_eur_v{version}_{dataset}.rg"
    conda: "../envs/ldsc.yaml"
    shell: "resources/ldsc/ldsc/ldsc.py --rg {input.mdd},{input.other} --ref-ld-chr {input.w_ld}/ --w-ld-chr {input.w_ld}/ --out {params.prefix}"

# return ids of GWAS datasets with SNP counts >= 2M and known sample sizes
# exclude a small number of studies that have non-sensical BETAs
def open_gwas_parse_dataset_ids(gwasinfo_txt, min_snps=2e6):
    import pandas as pd
    import os.path
    remove_ids = ['ebi-a-GCST007236', 'ebi-a-GCST002318', 'ebi-a-GCST006980', 'ebi-a-GCST007799', 'ebi-a-GCST007800']
    if(os.path.exists(gwasinfo_txt)):
        gwasinfo = pd.read_table(gwasinfo_txt)
        gwasinfo_pass = gwasinfo[(gwasinfo.nsnp >= min_snps) & (gwasinfo.sample_size.notna()) & (gwasinfo.id.isin(remove_ids) == False)]
        return(gwasinfo_pass.id.values)
    else:
        return([])
        
# extract table from LDSC log files
rule open_gwas_ldsc_log:
    input: "results/open_gwas/ldsc/{cohorts}_eur_v{version}_{dataset}.rg.log"
    output: "results/open_gwas/ldsc/{cohorts}_eur_v{version}_{dataset}.rg.txt"
    shell: "cat {input} | tail -n 5 | head -n 2 > {output}"
    
# list all IDs to fetch and run LDSC on, then generate report
rule open_gwas_rg_all_datasets:
    input: ldsc_full=expand("results/open_gwas/ldsc/full_eur_v{version}_{dataset}.rg.txt", version=analysis_version, dataset=open_gwas_parse_dataset_ids(rules.open_gwas_phewas_lookup.output.gwasinfo)), ldsc_noukbb=expand("results/open_gwas/ldsc/noUKBB_eur_v{version}_{dataset}.rg.txt", version=analysis_version, dataset=open_gwas_parse_dataset_ids(rules.open_gwas_phewas_lookup.output.gwasinfo)), gwasinfo="results/open_gwas/gwasinfo.txt"
    conda: "../envs/reports.yaml"
    output: full="docs/tables/ldsc_open_rg.txt", noukbb="docs/tables/ldsc_open_noUKBB.txt",  mr='docs/tables/ldsc_open_mr_candidates.txt'
    script: "../scripts/open/ldsc.R"
    
rule open_gwas:
    input: full="docs/tables/ldsc_open_rg.txt", noukbb="docs/tables/ldsc_open_noUKBB.txt",  mr='docs/tables/ldsc_open_mr_candidates.txt'
    conda: "../envs/reports.yaml"
    output: "docs/ldsc.md"
    script: "../docs/ldsc.Rmd"

