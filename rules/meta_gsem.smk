# Structured meta-analysis in GenomicSEM

# Grouping cohorts based on phenotypic measurement

# study types
# Clinical interview assessed (genotyped and sumstats)
# In-person interviews, clinic samples, review of medical records
meta_clin_geno=['antpo', 'bidi1', 'boma', 'cardm', 'cof3', 'col3', 'edi2', 'emcbp', 'formm', 'gens', 'gep3', 'grdg', 'grnd', 'gsk2', 'gsmse', 'gsrdf', 'gsrdg', 'gsrdi', 'gsrdp', 'i2b3', 'ihseu', 'mazdr', 'mmi2', 'mmo4', 'mrive', 'muen2', 'muspc', 'nes1', 'pfm2', 'qi3c', 'qi6c', 'qio2', 'rad3', 'rage', 'rai2', 'rau2', 'rde4', 'roc3', 'rot4', 'shp0', 'shpt', 'stm2', 'topmd', 'trail', 'twg2', 'yapeu']
meta_clin_sums=['GenScot', 'lgic2', 'tkda1']

# Health register/EHR
# Diagnostic record codes
meta_ehr_geno=['iruts']
meta_ehr_sums = ['iPSYCH', 'deCODE', 'HUNT', 'BioVU', 'PBK', 'SHARE', 'MoBa', 'MVP', 'GERA', 'DBDS', 'EXCEED', 'FinnGen', 'PREFECT', 'ESTBB']

# Questionnaire-derived diagnosis, completed by the participant
meta_quest_geno = ['prote']
meta_quest_sums = ['AGDS', 'ALSPAC', 'BASIC', 'STAGE', 'UKBB']

# Self-reported diagnosis
meta_selfrep_geno = []
meta_selfrep_sums = ['23andMe', 'Airwave']

meta_structured_groups_geno = {'Clin': meta_clin_geno, 'EHR': meta_ehr_geno, 'Quest': meta_quest_geno, 'SelfRep': meta_selfrep_geno}
meta_structured_groups_sums = {'Clin': meta_clin_sums, 'EHR': meta_ehr_sums, 'Quest': meta_quest_sums, 'SelfRep': meta_selfrep_sums}

# create reference info file linking to imputation panel
rule meta_gsem_refdir:
    output: "results/meta/gsem/reference_info"
    log: "logs/meta/gsem/reference_info.log"
    shell: "cd results/meta/gsem; impute_dirsub --refdir {config[refdir]} --reference_info --outname meta"
    
rule meta_gsem_sumstats:
    input: sumstats="results/sumstats/filtered/{cohort}.qc.gz"
    output: "results/meta/gsem/{cohort}.qc.gz"
    log: "logs/meta/{cohort}.log"
    shell: "ln -sv $(readlink -f {input.sumstats}) {output} > {log}"
    
if 'single_sumstats' in config:
    config_single_sumstats = config['single_sumstats']
else:
    config_single_sumstats = '~'
    
rule meta_gsem_single_sumstats:
    input: sumstats=expand("{single_sumstats}/daner_mdd_{{cohort}}_{{ancestries}}_{{qc}}.hg19.ch.fl.gz", single_sumstats=config_single_sumstats)
    output: "results/meta/gsem/daner_mdd_{cohort}_{ancestries}_{qc}.hg19.ch.fl.gz"
    shell: "ln -sv $(readlink -f {input.sumstats}) {output}"

def meta_gsem_cohorts(group, ancestry):
    cohorts = meta_structured_groups_sums[group]
    if(ancestry == 'eur'):
        cohorts_releases = [cohort for cohort in cohorts_eur if cohort[0] in cohorts]
    if(ancestry == 'eas'):
        cohorts_releases = [cohort for cohort in cohorts_eas if cohort[0] in cohorts]
    return(cohorts_releases)
    
    
def meta_gsem_cohorts_single(group, ancestry):
    cohorts = meta_structured_groups_geno[group]
    if(ancestry == 'eur'):
        cohorts_releases = [cohort for cohort in cohorts_geno_eur if cohort[0] in cohorts]
    if(ancestry == 'eas'):
        cohorts_releases = [cohort for cohort in cohorts_geno_eas if cohort[0] in cohorts]
    return(cohorts_releases)

# Ricopili results dataset list for eur structured 
rule meta_gsem_dataset_eur:
    input: lambda wildcards: expand("results/meta/gsem/daner_mdd_{cohort}_eur_{release}.hg19.ch.fl.gz", zip, cohort=[cohort[0] for cohort in meta_gsem_cohorts_single(wildcards.cohorts, 'eur')], release=[cohort[1] for cohort in meta_gsem_cohorts_single(wildcards.cohorts, 'eur')]), lambda wildcards: expand("results/meta/gsem/daner_mdd_{cohort}.eur.hg19.{release}.qc.gz", zip, cohort=[cohort[0] for cohort in meta_gsem_cohorts(wildcards.cohorts, 'eur')], release=[cohort[1] for cohort in meta_gsem_cohorts(wildcards.cohorts, 'eur')])
    output: "results/meta/gsem/dataset_{cohorts}_eur_v{version}"
    log: "logs/meta/gsem/dataset_{cohorts}_eur_v{version}.log"
    shell: "for daner in {input}; do echo $(basename $daner) >> {output}; done"
    
# Ricopili datasets files
rule meta_gsem_datasets_eur:
    input: expand("results/meta/gsem/dataset_{cohorts}_eur_v{version}", cohorts=meta_structured_groups_sums.keys(), version=analysis_version)
    
# Ricopili submission
rule meta_gsem_postimp:
    input: dataset="results/meta/gsem/dataset_{cohorts}_{ancestries}_v{version}", ref="results/meta/gsem/reference_info"
    params:
        popname=lambda wildcards: wildcards.ancestries.upper(),
        dataset=lambda wildcards, input: os.path.basename(input.dataset)
    output: touch("results/meta/gsem/{cohorts}_{ancestries}_v{version}.done")
    log: "logs/meta/gsem/pgc_mdd_meta_{cohorts}_{ancestries}_hg19_v{version}.postimp_navi.log"
    shell: "cd results/meta/gsem; postimp_navi --result {params.dataset} --popname {params.popname} --onlymeta --nolahunt --noldsc --out pgc_mdd_{wildcards.cohorts}_{wildcards.ancestries}_hg19_v{wildcards.version}"
	
rule install_gsem:
	output: "resources/ldsc/install_genomicsem.done"
	conda: "../envs/gsem.yaml"
	shell: """Rscript -e 'devtools::install_github("GenomicSEM/GenomicSEM", upgrade="never")' 2>&1 > {output}"""
    
rule install_ggmiami:
    output: "resources/ldsc/install_ggmiami.done"
    conda: "../envs/gsem.yaml"
    shell: """Rscript -e 'devtools::install_github("juliedwhite/miamiplot", upgrade="never")' 2>&1 > {output}"""
    
# sumstats to munge with Neff as sample size
# Neff as input sample size allows unbiased estimates of ldsc h2: https://www.medrxiv.org/content/10.1101/2021.09.22.21263909v1
rule meta_gsem_neff:
    input: "results/meta/gsem/distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz"
    output: "results/meta/gsem/neff/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.neff.txt"
    shell: """zcat {input} | awk '{{if(NR == 1) {{print $1, $2, $3, $4, $5, "MAF", $8, $9, $10, $11, "N"}} else {{print $1, $2, $3, $4, $5, $7, $8, $9, $10, $11, 2*$19}}}}' > {output}"""
    
rule meta_gsem_munge:
    input: sumstats="results/meta/gsem/neff/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.neff.txt", hm3="resources/ldsc/w_hm3.snplist", gsem="resources/ldsc/install_genomicsem.done"
    params: prefix="results/meta/gsem/munged/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}"
    output: "results/meta/gsem/munged/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.sumstats.gz"
    conda: "../envs/gsem.yaml"
    script: "../scripts/meta/gsem_munge.R"
    
	
# Create LDSC covstruct in GenomicSEM
rule meta_gsem_ldsc:
	input: sumstats=expand("results/meta/gsem/munged/pgc_mdd_{cohorts}_{{ancestries}}_hg19_v{version}.sumstats.gz", cohorts=meta_structured_groups_geno.keys(), version=analysis_version), samples=expand("results/meta/gsem/distribution/pgc_mdd_{cohorts}_{{ancestries}}_hg19_v{version}/basic.pgc_mdd_{cohorts}_eur_hg19_v{version}.num.xls", cohorts=meta_structured_groups_geno.keys(), version=analysis_version), w_ld_chr="resources/ldsc/{ancestries}_w_ld_chr/", gsem="resources/ldsc/install_genomicsem.done"
	params: cohorts=meta_structured_groups_geno.keys()
	output: covstruct="docs/objects/covstruct.{ancestries}.R", ldsc_table="docs/tables/meta_gsem_ldsc.{ancestries}.txt"
	conda: "../envs/gsem.yaml"
	script: "../scripts/meta/gsem_ldsc.R"
    
# Sumstats for usergwas
# Note: ref file downlaoded from https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v
rule meta_gsem_prepare_sumstats:
    input: sumstats=expand("results/meta/gsem/neff/pgc_mdd_{cohorts}_{{ancestries}}_hg19_v{{version}}.neff.txt", cohorts=meta_structured_groups_geno.keys()), ref="resources/gsem/reference.1000G.maf.0.005.txt"
    params: cohorts=meta_structured_groups_geno.keys()
    output: "results/meta/gsem/pgc_mdd_{ancestries}_v{version}.gsem.sumstats.gz"
    conda: "../envs/gsem.yaml"
    script: "../scripts/meta/gsem_sumstats.R"
    
# common factor GWAS, piece-wise
rule meta_gsem_commonfactor:
    input: sumstats="results/meta/gsem/pgc_mdd_{ancestries}_v{version}.gsem.sumstats.gz", covstruct="docs/objects/covstruct.{ancestries}.R"
    conda: "../envs/gsem.yaml"
    output: "results/meta/gsem/commonfactor/pgc_mdd_{ancestries}_v{version}.gsem.commonfactor.{CHR}-{START}-{STOP}.gz"
    script: "../scripts/meta/gsem_commonfactor.R"
    
def meta_gsem_1000G_regions(ref_file, k=1000):
    import os
    if(os.path.exists(ref_file)):
        import pandas as pd
        import numpy as np
        ref = pd.read_table(ref_file, sep='\s+')
        # remove chr23
        ref_auto = ref[ref['CHR'] != 23].copy() 
        n_snps = ref_auto.shape[0]
        # index to group SNPs into regions
        region_index = np.arange(n_snps) // (n_snps // k)
        ref_auto['region'] = region_index[0:n_snps]
        # get start and stop coordinates of each region
        regions = ref_auto.groupby(['CHR', 'region']).agg(start=pd.NamedAgg(column="BP", aggfunc="min"), stop=pd.NamedAgg(column="BP", aggfunc="max"))
        # turn grouping indices into columns
        regions.reset_index(inplace=True)
        # concatenate list of CHR-START-STOP regions
        regions_str = regions['CHR'].map(str) + '-' + regions['start'].map(str) + '-' + regions['stop'].map(str)
    else:
        regions_str = ''
    return(regions_str)

# merge common factor GWAS
rule meta_gsem_commonfactor_merge:
    input: sumstats=expand("results/meta/gsem/commonfactor/pgc_mdd_{{ancestries}}_v{{version}}.gsem.commonfactor.{region}.gz", region=meta_gsem_1000G_regions(ref_file="resources/gsem/reference.1000G.maf.0.005.txt"))
    output: "results/meta/gsem/pgc_mdd_{ancestries}_v{version}.gsem.commonfactor.gz"
    conda: "../envs/gsem.yaml"
    script: "../scripts/meta/gsem_commonfactor_merge.R"
    
# reformat as daner
rule meta_gsem_commonfactor_daner:
    input: gsem="results/meta/gsem/pgc_mdd_{ancestries}_v{version}.gsem.commonfactor.gz", daner="results/distribution/daner_pgc_mdd_full_{ancestries}_hg19_v{version}.rp.gz"
    output: "results/meta/gsem/daner_pgc_mdd_full_{ancestries}_v{version}.gsem.commonfactor.gz"
    conda: "../envs/gsem.yaml"
    script: "../scripts/meta/gsem_daner.R"
    
# Notebook
rule meta_gsem:
    input: covstruct_eur="docs/objects/covstruct.eur.R", ldsc_table_eur="docs/tables/meta_gsem_ldsc.eur.txt", commonfactor_eur=expand("results/meta/gsem/pgc_mdd_eur_v{version}.gsem.commonfactor.gz", version=analysis_version), notebook="docs/gsem.Rmd", gsem="resources/ldsc/install_genomicsem.done"
    output: "docs/gsem.md"
    conda: "../envs/gsem.yaml"
    script: "../docs/gsem.Rmd"