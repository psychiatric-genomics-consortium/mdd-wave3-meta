
rule install_polyfun: 
    output: directory("resources/finemapping/polyfun") 
    conda: "../envs/finemapping.yaml" 
    shell: "git clone git@github.com:omerwe/polyfun.git {output}" 

rule get_baseline
  input: HTTP.remote("//data.broadinstitute.org/alkesgroup/LDSCORE/baselineLF_v2.2.UKB.polyfun.tar.gz")
  output: "resources/finemapping/baselineLF2.2.UKB.polyfun.tar.gz"
  shell: "cp {input} {output}; tar -xzvf {output}"

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

rule define_n
    input: "results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz"
    output: "results/finemapping/n_{cohorts}_{ancestries}_hg19_v{version}.gz"
    shell:
        "declare -i cases controls n
        cases=$(gunzip -c {input} | head -1 | awk '{print $6}' | sed 's/FRQ_A_//g')
        controls=$(gunzip -c {input} | head -1 | awk '{print $7}' | sed 's/FRQ_U_//g')
        n=$cases+$controls
        echo $n > {output}"

rule create_and_run_jobs:
    input:
        snpvar=expand("results/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}_chr{chr}.snpvar_ridge_constrained.gz", chr=range(1, 22)),
        n="results/finemapping/n_{cohorts}_{ancestries}_hg19_v{version}.gz
    output: "results/finemapping/results_{cohorts}_{ancestries}_hg19_v{version}"
    script: "../scripts/finemapping/create.sh {input.snpvar} {input.n} {output}"
