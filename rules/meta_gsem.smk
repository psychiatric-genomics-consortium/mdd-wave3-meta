# Structured meta-analysis in GenomicSEM

# study types
# Clinically assessed
meta_clin=['MDD49', 'tkda1', 'lgic2']

# Health register/EHR
meta_ehr = ['iPSYCH', 'HUNT', 'BioVU', 'PBK', 'SHARE', 'MoBa', 'MVP', 'GERA', 'DBDS', 'EXCEED', 'FinnGen', 'PREFECT', 'ESTBB']

# mixed
meta_mixed = ['GenScot', 'UKBB', 'deCODE']

# questionnaire
meta_quest = ['AGDS', 'ALSPAC', 'BASIC', 'STAGE']

# self-declared
meta_self = ['23andMe', 'Airwave']

meta_structured_groups = {'clin': meta_clin, 'ehr': meta_ehr, 'mixed': meta_mixed, 'quest': meta_quest, 'self': meta_self}

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
    
rule meta_gsem_datasets_eur:
    input: expand("results/meta/gsem/dataset_{cohorts}_eur_v{version}", cohorts=meta_structured_groups.keys(), version=analysis_version)