rule install_polyfun: 
    output: directory("resources/finemapping/polyfun") 
    conda: "../envs/finemapping.yaml"
    log: "logs/finemapping/install_polyfun.log"
    shell: "git clone git@github.com:omerwe/polyfun.git {output}" 

rule install_susie:
    conda: "../envs/finemapping.yaml"
    log: "logs/finemapping/install_susie.log"
    script: "../scripts/finemapping/install_susie.R"

rule get_baseline:
    input: HTTP.remote("data.broadinstitute.org/alkesgroup/LDSCORE/baselineLF_v2.2.UKB.polyfun.tar.gz", keep_local=False)
    log: "logs/finemapping/get_baseline.log"
    conda: "../envs/finemapping.yaml" 
    output: "resources/finemapping/baselineLF_v2.2.UKB.polyfun.tar.gz"
    shell: "tar -xzvf {input}; mv baselineLF2.2.UKB resources/finemapping/"

rule format_sumstat:
    input: "results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz"
    output: "results/finemapping/for_munge_{cohorts}_{ancestries}_hg19_v{version}.gz"
    log: "logs/finemapping/for_munge_{cohorts}_{ancestries}_hg19_v{version}.log"
    conda: "../envs/finemapping.yaml"
    shell: "scripts/finemapping/format.bash {input} {output}"

rule munge_sumstat:
    input: sumstats="results/finemapping/for_munge_{cohorts}_{ancestries}_hg19_v{version}.gz",
        polyfun="resources/finemapping/polyfun"
    output: "results/finemapping/munged_{cohorts}_{ancestries}_hg19_v{version}.parquet"
    log: "logs/finemapping/munged_{cohorts}_{ancestries}_hg19_v{version}.log"
    conda: "../envs/finemapping.yaml"
    shell:
        "python3 {input.polyfun}/munge_polyfun_sumstats.py \
        --sumstats {input.sumstats} \
        --out {output} \
        --min-info 0.6 \
        --min-maf 0.001"

rule l2reg_sldsc:
    input: sumstats="results/finemapping/munged_{cohorts}_{ancestries}_hg19_v{version}.parquet",
        polyfun="resources/finemapping/polyfun",
	baseline="logs/finemapping/get_baseline.log"
    params: out="results/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}_chr"
    log: "logs/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}.log"
    conda: "../envs/finemapping.yaml" 
    shell:
        "python3 {input.polyfun}/polyfun.py \
        --compute-h2-L2 \
        --no-partitions \
        --output-prefix {params.out} \
        --sumstats {input.sumstats} \
        --ref-ld-chr resources/finemapping/baselineLF2.2.UKB/baselineLF2.2.UKB. \
        --w-ld-chr resources/finemapping/baselineLF2.2.UKB/weights.UKB. \
        --allow-missing"

rule define_n:
    input: "results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz"
    output: "results/finemapping/n_{cohorts}_{ancestries}_hg19_v{version}.out"
    log: "logs/finemapping/n_{cohorts}_{ancestries}_hg19_v{version}.log"
    conda: "../envs/finemapping.yaml" 
    shell: "scripts/finemapping/define_n.bash {input} {output}"

rule create_and_run_jobs:
    input: snpvar="logs/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}.log",
        n="results/finemapping/n_{cohorts}_{ancestries}_hg19_v{version}.out",
	susie="logs/finemapping/install_susie.log"
    params: snpvar="results/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}_chr"
    output: "results/finemapping/results_{cohorts}_{ancestries}_hg19_v{version}"
    log: "logs/finemapping/results_{cohorts}_{ancestries}_hg19_v{version}.log"
    conda: "../envs/finemapping.yaml"
    shell: "scripts/finemapping/create.sh {params.snpvar} {input.n} {output}"
