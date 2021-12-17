###################################
# Conduct meta-analysis in Ricopili
###################################

###
### Cohorts sets
###
    
# cohort sets for analysts
cohorts_analyst = ["full", "noUKBB", "noALSPAC"]

# cohort sets for public
cohorts_public = ["no23andMe"]

###
### Split and X
###

# split sumstats files into autosome and X chromosome
rule meta_chr22:
    input: "results/sumstats/filtered/{cohort}.gz"
    output: "results/sumstats/autosome/{cohort}.gz"
    shell: "zcat {input} | awk '{{if(NR == 1 || $1 != 23) {{print $0}}}}' | gzip -c > {output}"
    
rule meta_chr23:
    input: "results/sumstats/filtered/{cohort}.gz"
    output: "results/sumstats/X/{cohort}.gz"
    shell: "zcat {input} | awk '{{if(NR == 1 || $1 == 23) {{print $0}}}}' | gzip -c > {output}"

###
### Link sumstats inside meta analysis directory
###

# link autosome sumstats files into meta-analysis directory, but also run
# LDSC rg with MDD2
rule meta:
	input: sumstats="results/sumstats/autosome/{cohort}.gz", rg="results/sumstats/rg_mdd/{cohort}.log"
	output: "results/meta/{cohort}.gz"
	log: "logs/meta/{cohort}.log"
	shell: "ln -sv $(readlink -f {input.sumstats}) {output} > {log}"

# link chrX sumstats files into meta-analysis directory
rule metax:
    input: sumstats="results/sumstats/X/{cohort}.gz"
    output: "results/meta/X/{cohort}.gz"
    log: "logs/meta/{cohort}.log"
    shell: "ln -sv $(readlink -f {input.sumstats}) {output} > {log}"

###
### List sumstats in Ricopili dataset files
###

# Ricopili results dataset list for eur ancestries
rule dataset_eur:
	input: expand("results/meta/daner_mdd_{cohort}.eur.hg19.{release}.qc.gz", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur])
	output: "results/meta/dataset_full_eur_v{analysis}"
	log: "logs/meta/dataset_full_eur_v{analysis}.log"
	shell: "for daner in {input}; do echo $(basename $daner) >> {output}; done"
	
# Ricopili results dataset list for eas ancestries
rule dataset_eas:
	input: expand("results/meta/daner_mdd_{cohort}.eas.hg19.{release}.qc.gz", zip, cohort=[cohort[0] for cohort in cohorts_eas], release=[cohort[1] for cohort in cohorts_eas])
	output: "results/meta/dataset_full_eas_v{analysis}"
	log: "logs/meta/dataset_full_eas_v{analysis}.log"
	shell: "for daner in {input}; do echo $(basename $daner) >> {output}; done"
	
# Dataset list that exclude a particular cohort
rule dataset_noCOHORT_eur:
	input: "results/meta/dataset_full_eur_v{analysis}"
	log: "logs/meta/dataset_no{cohort}_eur_v{analysis}"
	output: "results/meta/dataset_no{cohort}_eur_v{analysis}"
	shell: "cat {input} | grep --invert daner_mdd_{wildcards.cohort} > {output}"
    
# Ricopili results dataset list for eur ancestries chrX
rule dataset_eur_X:
    input: expand("results/meta/X/daner_mdd_{cohort}.eur.hg19.{release}.qc.gz", zip, cohort=[cohort[0] for cohort in cohorts_eur], release=[cohort[1] for cohort in cohorts_eur])
    output: "results/meta/X/dataset_full_eur_v{analysis}"
    log: "logs/meta/X/dataset_full_eur_v{analysis}.log"
    shell: """for daner in {input}; do 
    headn=$(zcat $daner | awk '$1 == 23' | head | wc -l)
    if ((headn > 0)); then
        echo $(basename $daner) >> {output}; 
    fi;
    done
    """
    
# Dataset list that exclude a particular X cohort
rule dataset_noCOHORT_X_eur:
    input: "results/meta/X/dataset_full_eur_v{analysis}"
    log: "logs/meta/X/dataset_no{cohort}_eur_v{analysis}"
    output: "results/meta/X/dataset_no{cohort}_eur_v{analysis}"
    shell: "cat {input} | grep --invert daner_mdd_{wildcards.cohort} > {output}"
    
rule dataset_eas_X:
    input: expand("results/meta/X/daner_mdd_{cohort}.eas.hg19.{release}.qc.gz", zip, cohort=[cohort[0] for cohort in cohorts_eas], release=[cohort[1] for cohort in cohorts_eas])
    output: "results/meta/X/dataset_full_eas_v{analysis}"
    log: "logs/meta/X/dataset_full_eas_v{analysis}.log"
    shell: """for daner in {input}; do 
    headn=$(zcat $daner | awk '$1 == 23' | head | wc -l)
    if ((headn > 0)); then
        echo $(basename $daner) >> {output}; 
    fi;
    done
    """

###
### Run Ricopili
###

# Ricopili submission
rule postimp:
	input: dataset="results/meta/dataset_{cohorts}_{ancestries}_v{version}", ref="results/meta/reference_info"
	params:
		popname=lambda wildcards: wildcards.ancestries.upper(),
		dataset=lambda wildcards, input: os.path.basename(input.dataset)
	output: touch("results/meta/{cohorts}_{ancestries}_v{version}.done")
	log: "logs/meta/pgc_mdd_meta_{cohorts}_{ancestries}_hg19_v{version}.postimp_navi.log"
	shell: "cd results/meta; postimp_navi --result {params.dataset} --popname {params.popname} --nolahunt --noldsc --no_neff_filter --out pgc_mdd_{wildcards.cohorts}_{wildcards.ancestries}_hg19_v{wildcards.version}"

rule postimp_eur:
	input: expand("results/meta/full_eur_v{version}.done", version=analysis_version_eur)
	
rule postimp_eas:
	input: expand("results/meta/full_eas_v{version}.done", version=analysis_version_eas)
	
# Ricopili submission chrX
rule postimpX:
    input: dataset="results/meta/X/dataset_{cohorts}_{ancestries}_v{version}", ref="results/meta/X/reference_info"
    params:
        popname=lambda wildcards: wildcards.ancestries.upper(),
        dataset=lambda wildcards, input: os.path.basename(input.dataset)
    output: touch("results/meta/X/{cohorts}_{ancestries}_v{version}.done")
    log: "logs/meta/X/pgc_mdd_meta_{cohorts}_{ancestries}_hg19_v{version}.postimp_navi.log"
    shell: "cd results/meta/X; postimp_navi --result {params.dataset} --popname {params.popname} --nolahunt --noldsc --no_neff_filter --out pgc_mdd_{wildcards.cohorts}_{wildcards.ancestries}_hg19_v{wildcards.version}"
    
    
rule postimp_eur_X:
    input: expand("results/meta/X/full_eur_v{version}.done", version=analysis_version)
        
    
###
### Clump results with custom QC params
###
    
# Clump results based on workflow params
rule postimp_reclump:
	input: "results/meta/distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz"
	params:
		refdir=config['refdir'],
		p1=meta_qc_params['clu_p1'],
		p2=meta_qc_params['clu_p2'],
		r2=meta_qc_params['clu_r2'],
		frq=meta_qc_params['clu_maf'],
		info=meta_qc_params['clu_info'],
		window=meta_qc_params['clu_kb'],
		popname=lambda wildcards: wildcards.ancestries.upper(),
		pfile="distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz",
		outname="pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.reclump"
	output: touch("results/meta/{cohorts}_{ancestries}_v{version}.reclump.done")
	shell: "cd results/meta; clump_nav3 --pfile {params.pfile} --refdir {params.refdir}/pop_{params.popname} --clu_p1 {params.p1} --clu_p2 {params.p2} --clu_r2 {params.r2} --clu_window {params.window} --popname {params.popname} --outname {params.outname} --debug --serial --sepa 16"
	
    # clump for top 10k results
rule postimp_reclump_top10k:
    input: "results/meta/distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.neff.gz"
    params:
        refdir=config['refdir'],
        p1=meta_qc_params['clu10k_p1'],
        p2=meta_qc_params['clu10k_p2'],
        r2=meta_qc_params['clu10k_r2'],
        frq=meta_qc_params['clu10k_maf'],
        info=meta_qc_params['clu10k_info'],
        window=meta_qc_params['clu10k_kb'],
        popname=lambda wildcards: wildcards.ancestries.upper(),
        pfile="distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.neff.gz",
        outname="pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.top10k"
    output: touch("results/meta/{cohorts}_{ancestries}_v{version}.top10k.done")
    shell: """
    cd results/meta; 
    clump_nav3 --pfile {params.pfile} \
    --refdir {params.refdir}/pop_{params.popname} \
    --clu_p1 {params.p1} \
    --clu_p2 {params.p2} \
    --clu_r2 {params.r2} \
    --clu_window {params.window} \
    --hq_f {params.frq} \
    --hq_i {params.info} \
    --popname {params.popname} \
    --outname {params.outname} \
    --debug --serial --sepa 16
    """

