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


# Ricopili results dataset list for eur structured 
rule meta_gsem_dataset_eur:
    input: expand("results/meta/daner_mdd_{cohort}.eur.hg19.{release}.qc.gz", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur])
    output: "results/meta/dataset_full_eur_v{analysis}"
    log: "logs/meta/dataset_full_eur_v{analysis}.log"
    shell: "for daner in {input}; do echo $(basename $daner) >> {output}; done"