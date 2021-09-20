# Structured meta-analysis in GenomicSEM

# study types
# Clinical interview assessed (genotyped and sumstats)
meta_clin_geno=['gep3', 'grdg', 'grnd', 'gsk2', 'gsmse', 'gsrdf', 'gsrdg', 'gsrdi', 'gsrdp', 'i2b3', 'ihseu', 'jjp2', 'mazdr', 'mmi2', 'mmo4', 'mrive', 'muen2', 'muspc', 'nes1', 'pfm2', 'qi3c', 'qi6c', 'qio2', 'rad3', 'rage', 'rai2', 'rau2', 'rde4', 'roc3', 'rot4', 'shp0', 'shpt', 'stm2', 'topmd', 'trail', 'twg2', 'yapeu']
meta_clin_sums=['GenScot', 'lgic2', 'tkda1']

# Health register/EHR
meta_ehr_geno=['iruts']
meta_ehr_sums = ['iPSYCH', 'deCODE', 'HUNT', 'BioVU', 'PBK', 'SHARE', 'MoBa', 'MVP', 'GERA', 'DBDS', 'EXCEED', 'FinnGen', 'PREFECT', 'ESTBB']

# questionnaire derived diagnosis
meta_quest_geno = ['prote']
meta_quest_sums = ['AGDS', 'ALSPAC', 'BASIC', 'STAGE', 'UKBB']

# self-reported diagnosis
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
    
rule meta_gsem_single_sumstats:
    input: sumstats=expand("{single_sumstats}/daner_mdd_{{cohort}}_{{ancestries}}_{{qc}}.hg19.ch.fl.gz", single_sumstats=config['single_sumstats'])
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
    shell: "cd results/meta/gsem; postimp_navi --result {params.dataset} --popname {params.popname} --onlymeta --nolahunt --out pgc_mdd_{wildcards.cohorts}_{wildcards.ancestries}_hg19_v{wildcards.version}"
	
rule install_gsem:
	output: "resources/ldsc/install_genomicsem.done"
	conda: "../envs/gsem.yaml"
	shell: """Rscript -e 'devtools::install_github("GenomicSEM/GenomicSEM", upgrade="never")' 2>&1 > {output}"""
	
# Create LDSC covstruct in GenomicSEM
rule meta_gsem_ldsc:
	input: sumstats=expand("results/meta/gsem/distribution/pgc_mdd_{cohorts}_{{ancestries}}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{{ancestries}}_hg19_v{version}.gz.ldsc.sumstats.gz", cohorts=meta_structured_groups_geno.keys(), version=analysis_version), samples=expand("results/meta/gsem/distribution/pgc_mdd_{cohorts}_{{ancestries}}_hg19_v{version}/basic.pgc_mdd_{cohorts}_eur_hg19_v{version}.num.xls", cohorts=meta_structured_groups_geno.keys(), version=analysis_version), w_ld_chr="resources/ldsc/{ancestries}_w_ld_chr/", gsem="resources/ldsc/install_genomicsem.done"
	params: cohorts=meta_structured_groups_geno.keys()
	output: covstruct="docs/objects/covstruct.{ancestries}.R", ldsc_table="docs/tables/meta_gsem_ldsc.{ancestries}.txt"
	conda: "../envs/gsem.yaml"
	script: "../scripts/meta/gsem_ldsc.R"
    
# Notebook
rule meta_gsem:
    input: covstruct_eur="docs/objects/covstruct.eur.R", ldsc_table_eur="docs/tables/meta_gsem_ldsc.eur.txt", notebook="docs/gsem.Rmd", gsem="resources/ldsc/install_genomicsem.done"
    output: "docs/gsem.md"
    conda: "../envs/gsem.yaml"
    script: "../docs/gsem.Rmd"