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
    
# QC based on INFO score and Neff
rule postimp_rp_neff:
    input: "results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.rp.gz"
    log: "logs/meta/distribution/{analysis}.neff.log"
    params: qc_neff=0.8
    conda: "../envs/meta.yaml"
    output: "results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.neff.gz"
    script: "../scripts/meta/rp_neff.R"

# inputs for postimp_rp
rule postimp_rp_all:
    input: expand("results/meta/distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.neff.gz", cohorts=cohorts_analyst + cohorts_public, ancestries=['eur'], version=analysis_version)