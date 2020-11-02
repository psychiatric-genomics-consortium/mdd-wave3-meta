
rule format_sumstat:
    input: "results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz"
    output: "results/finemapping/for_munge_{cohorts}_{ancestries}_hg19_v{version}.gz
    script: "../scripts/finemapping/format.bash {input} {output}"

rule munge_sumstat:
    input: "results/finemapping/for_munge_{cohorts}_{ancestries}_hg19_v{version}.gz
    output: "results/finemapping/munged_{cohorts}_{ancestries}_hg19_v{version}.parquet"
    shell:
        "python3 resources/finemapping/polyfun/munge_polyfun_sumstats.py \
        --sumstats {input} \
        --out {output} \
        --min-info 0.6 \
        --min-maf 0.001"

rule l2reg_sldsc:
    input: "results/finemapping/munged_{cohorts}_{ancestries}_hg19_v{version}.parquet"
    output: "results/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}_chr"
    shell:
        "python3 resources/finemapping/polyfun/munge_polyfun_sumstats.py \
        --compute-h2-L2 \
        --no-partitions \
        --output-prefix {output} \
        --sumstats {input} \
        --ref-ld-chr resources/finemapping/baselineLF2.2.UKB/baselineLF2.2.UKB. \
        --w-ld-chr resources/finemapping/baselineLF2.2.UKB/weights.UKB. \
        --allow-missing"

### NOTE TO FUTURE JONI - AUTOMATE N BELOW

rule create_and_run_jobs:
    input:
        snpvar=expand("results/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}_chr{chr}.snpvar_ridge_constrained.gz", chr=range(1, 22)),
        n=12345
    output: "results/finemapping/results_{cohorts}_{ancestries}_hg19_v{version}"
    script: "../scripts/finemapping/create.sh {input.snpvar} {input.n} {output}"

