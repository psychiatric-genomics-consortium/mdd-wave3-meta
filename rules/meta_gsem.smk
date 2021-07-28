# Structured meta-analysis in GenomicSEM

# study types
# Clinically assessed
meta_clin=['MDD49', 'GenScot', 'lgic2', 'tkda1']

# Health register/EHR
meta_ehr = ['iPSYCH', 'deCODE', 'HUNT', 'BioVU', 'PBK', 'SHARE', 'MoBa', 'MVP', 'GERA', 'DBDS', 'EXCEED', 'FinnGen', 'PREFECT', 'ESTBB']

# questionnaire
meta_quest = ['AGDS', 'ALSPAC', 'BASIC', 'STAGE', 'UKBB']

# self-declared
meta_self = ['23andMe', 'Airwave']

meta_structured_groups = {'clin': meta_clin, 'ehr': meta_ehr, 'quest': meta_quest, 'self': meta_self}

# create reference info file linking to imputation panel
rule meta_gsem_refdir:
    output: "results/meta/gsem/reference_info"
    log: "logs/meta/gsem/reference_info.log"
    shell: "cd results/meta/gsem; impute_dirsub --refdir {config[refdir]} --reference_info --outname meta"
    
rule meta_gsem_sumstats:
    input: sumstats="results/sumstats/filtered/{cohort}.gz"
    output: "results/meta/gsem/{cohort}.gz"
    log: "logs/meta/{cohort}.log"
    shell: "ln -sv $(readlink -f {input.sumstats}) {output} > {log}"

def meta_gsem_cohorts(group, ancestry):
    cohorts = meta_structured_groups[group]
    if(ancestry == 'eur'):
        cohorts_releases = [cohort for cohort in cohorts_eur if cohort[0] in cohorts]
    if(ancestry == 'eas'):
        cohorts_releases = [cohort for cohort in cohorts_eas if cohort[0] in cohorts]
    return(cohorts_releases)

# Ricopili results dataset list for eur structured 
rule meta_gsem_dataset_eur:
    input: lambda wildcards: expand("results/meta/gsem/daner_mdd_{cohort}.eur.hg19.{release}.qc.gz", zip, cohort=[cohort[0] for cohort in meta_gsem_cohorts(wildcards.cohorts, 'eur')], release=[cohort[1] for cohort in meta_gsem_cohorts(wildcards.cohorts, 'eur')])
    output: "results/meta/gsem/dataset_{cohorts}_eur_v{version}"
    log: "logs/meta/gsem/dataset_{cohorts}_eur_v{version}.log"
    shell: "for daner in {input}; do echo $(basename $daner) >> {output}; done"
    
# Ricopili datasets files
rule meta_gsem_datasets_eur:
    input: expand("results/meta/gsem/dataset_{cohorts}_eur_v{version}", cohorts=meta_structured_groups.keys(), version=analysis_version)
    
# Ricopili submission
rule meta_gsem_postimp:
    input: dataset="results/meta/gsem/dataset_{cohorts}_{ancestries}_v{version}", ref="results/meta/gsem/reference_info"
    params:
        popname=lambda wildcards: wildcards.ancestries.upper(),
        dataset=lambda wildcards, input: os.path.basename(input.dataset)
    output: touch("results/meta/gsem/{cohorts}_{ancestries}_v{version}.done")
    log: "logs/meta/gsem/pgc_mdd_meta_{cohorts}_{ancestries}_hg19_v{version}.postimp_navi.log"
    shell: "cd results/meta/gsem; postimp_navi --result {params.dataset} --popname {params.popname} --onlymeta --nolahunt --no_neff_filter --out pgc_mdd_{wildcards.cohorts}_{wildcards.ancestries}_hg19_v{wildcards.version}"
	
rule install_gsem:
	output: "resources/ldsc/install_genomicsem.done"
	conda: "../envs/gsem.yaml"
	shell: """Rscript -e 'devtools::install_github("GenomicSEM/GenomicSEM", upgrade="never")' 2>&1 > {output}"""
	
# Create LDSC covstruct in GenomicSEM
rule meta_gsem_ldsc:
	input: sumstats=expand("results/meta/gsem/distribution/pgc_mdd_{cohorts}_{{ancestries}}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{{ancestries}}_hg19_v{version}.gz.ldsc.sumstats.gz", cohorts=meta_structured_groups.keys(), version=analysis_version), samples=expand("results/meta/gsem/distribution/pgc_mdd_{cohorts}_{{ancestries}}_hg19_v{version}/basic.pgc_mdd_{cohorts}_eur_hg19_v{version}.num.xls", cohorts=meta_structured_groups.keys(), version=analysis_version), w_ld_chr="resources/ldsc/{ancestries}_w_ld_chr/", gsem="resources/ldsc/install_genomicsem.done"
	params: cohorts=meta_structured_groups.keys()
	output: covstruct="docs/objects/covstruct.{ancestries}.R", ldsc_table="docs/tables/meta_gsem_ldsc.{ancestries}.txt"
	conda: "../envs/gsem.yaml"
	script: "../scripts/meta/gsem_ldsc.R"
    
# Notebook
rule meta_gsem:
    input: covstruct_eur="docs/objects/covstruct.eur.R", ldsc_table_eur="docs/tables/meta_gsem_ldsc.eur.txt", notebook="docs/gsem.Rmd", gsem="resources/ldsc/install_genomicsem.done"
    output: "docs/gsem.md"
    conda: "../envs/gsem.yaml"
    script: "../docs/gsem.Rmd"