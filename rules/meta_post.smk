############################################
## Post processing of meta-analysed sumstats
############################################

###
### Post processing
###

# check Ricopili output for complete rows and duplicates
rule postimp_rp:
    input: autosome="results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.gz", xsome="results/meta/X/report_pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.gz"
    log: "logs/meta/distribution/{analysis}.rp.log"
    conda: "../envs/meta.yaml"
    output: "results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.rp.gz"
    script: "../scripts/meta/rp.R"
    
# QC based on Neff
rule postimp_rp_neff:
    input: "results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.rp.gz"
    log: "logs/meta/distribution/{analysis}.neff.log"
    params: qc_neff=0.8
    conda: "../envs/meta.yaml"
    output: "results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.neff.gz"
    script: "../scripts/meta/rp_neff.R"

# Merge auto- and allosome clumped results
rule postimp_clumped:
    input: auto="results/meta/distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz.p4.clump.areator.sorted.1mhc", allo="results/meta/X/report_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz.p4.clump.areator.sorted.1mhc"
    output: "results/meta/distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz.p4.clump.areator.sorted.1mhc.txt"
    conda: "../envs/meta.yaml"
    script: "../scripts/meta/rp_clump.R"

# inputs for postimp_rp
rule postimp_rp_cohorts:
    input: neff=expand("results/meta/distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.neff.gz", cohorts=cohorts_analyst + cohorts_public, ancestries=['eur'], version=analysis_version), clump=expand("results/meta/distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz.p4.clump.areator.sorted.1mhc.txt", cohorts=['full'], ancestries=['eur'], version=analysis_version)
    
# add clumped results file to repo
rule postimp_rp_clumped:
    input: expand("results/meta/distribution/pgc_mdd_full_{{ancestries}}_hg19_v{version}/daner_pgc_mdd_full_{{ancestries}}_hg19_v{version}.gz.p4.clump.areator.sorted.1mhc.txt", version=analysis_version)
    output: "docs/tables/meta_snps_full_{ancestries}.clump.txt"
    shell: "cp {input} {output}"
    
# filter top10k clumped results
rule postimp_rp_top10k:
    input: auto="results/meta/pgc_mdd_full_{ancestries}_hg19_v{version}.top10k.clumped.xmhc.gz", allo="results/meta/X/pgc_mdd_full_{ancestries}_hg19_v{version}.neff.top10k.clumped", daner="results/meta/distribution/pgc_mdd_full_{ancestries}_hg19_v{version}/daner_pgc_mdd_full_{ancestries}_hg19_v{version}.rp.gz", cojo="docs/tables/meta_snps_full_{ancestries}.cojo.txt", clump="docs/tables/meta_snps_full_{ancestries}.clump.txt" 
    output: "results/meta/distribution/pgc_mdd_top10k_{ancestries}_hg19_v{version}/daner_pgc_mdd_top10k_{ancestries}_hg19_v{version}.neff.gz"
    conda: "../envs/meta.yaml"
    script: "../scripts/meta/rp_top10k.R"
    
rule postimp_rp_clumped_txt:
    input: "docs/tables/meta_snps_full_eur.clump.txt"