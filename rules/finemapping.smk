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
    shell: "tar -xzvf {input}; mv baselineLF2.2.UKB resources/finemapping/"

rule format_sumstat:
    input: "results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.rp.gz"
    output: "results/finemapping/for_munge_{cohorts}_{ancestries}_hg19_v{version}.rp.gz"
    log: "logs/finemapping/for_munge_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    conda: "../envs/finemapping.yaml"
    shell: "bash scripts/finemapping/format.bash {input} {output}"

rule munge_sumstat:
    input: sumstats="results/finemapping/for_munge_{cohorts}_{ancestries}_hg19_v{version}.rp.gz",
        polyfun="resources/finemapping/polyfun"
    output: "results/finemapping/munged_{cohorts}_{ancestries}_hg19_v{version}.rp.parquet"
    log: "logs/finemapping/munged_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    conda: "../envs/finemapping.yaml"
    shell:
        "python3 {input.polyfun}/munge_polyfun_sumstats.py \
        --sumstats {input.sumstats} \
        --out {output} \
        --min-info 0.6 \
        --min-maf 0.001 "

rule l2reg_sldsc:
    input: sumstats="results/finemapping/munged_{cohorts}_{ancestries}_hg19_v{version}.rp.parquet",
        polyfun="resources/finemapping/polyfun",
	baseline="logs/finemapping/get_baseline.log"
    params: out="results/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}.rp_chr"
    log: "logs/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    conda: "../envs/finemapping.yaml" 
    shell:
        "export PYTHONNOUSERSITE=True && python3 {input.polyfun}/polyfun.py \
        --compute-h2-L2 \
        --no-partitions \
        --output-prefix {params.out} \
        --sumstats {input.sumstats} \
        --ref-ld-chr resources/finemapping/baselineLF2.2.UKB/baselineLF2.2.UKB. \
        --w-ld-chr resources/finemapping/baselineLF2.2.UKB/weights.UKB. \
        --allow-missing"

rule define_n:
    input: "results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.rp.gz"
    output: "results/finemapping/n_{cohorts}_{ancestries}_hg19_v{version}.rp.out"
    log: "logs/finemapping/n_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    conda: "../envs/finemapping.yaml" 
    shell: "bash scripts/finemapping/define_n.bash {input} {output}"

rule create_genomewide_finemapping_jobs:
    input: snpvar="logs/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}.rp.log",
        n="results/finemapping/n_{cohorts}_{ancestries}_hg19_v{version}.rp.out",
	susie="logs/finemapping/install_susie.log"
    params: snpvar="results/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}.rp_chr",
        outprefix="results/finemapping/results_{cohorts}_{ancestries}_hg19_v{version}.rp"
    output: expand("results/finemapping/results_{{cohorts}}_{{ancestries}}_hg19_v{{version}}.rp_jobs_chr{chr}.sh", chr=range(1,22,1)) 
    log: "logs/finemapping/create_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    conda: "../envs/finemapping.yaml"
    shell: "export PYTHONNOUSERSITE=True && scripts/finemapping/create.sh {params.snpvar} {input.n} {params.outprefix}"

# RULE HASHED - COVERED BY CREATE.SH ABOVE
#
# rule run_finemapping_jobs_main:
#    input: expand("results/finemapping/results_{{cohorts}}_{{ancestries}}_hg19_v{{version}}.rp_jobs_chr{chr}.sh", chr=range(1,22,1))
#    output: expand("results/finemapping/results_{{cohorts}}_{{ancestries}}_hg19_v{{version}}.rp.chr{chr}*.gz", chr=range(1, 22, 1))
#    log: expand("logs/finemapping/run_{{cohorts}}_{{ancestries}}_hg19_v{{version}}.rp_jobs_chr{chr}.log", chr=range(1, 22, 1))
#    conda: "../envs/finemapping.yaml"
#    shell: "sh {input}"

rule merge_finemapping_jobs:
    input: expand("results/finemapping/results_{{cohorts}}_{{ancestries}}_hg19_v{{version}}.rp.chr{chr}*.gz", chr=range(1, 22, 1))
    output: "results/finemapping/merged_results_{cohorts}_{ancestries}_hg19_v{version}.rp_WG_susie_1.gz"
    log: "logs/finemapping/merged_results_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    conda: "../envs/finemapping.yaml" 
    shell: "cat ${input} >> ${output} && gunzip -c ${output} | sort -k11,11gr | head | column -t"

rule get_genomewide_loci:
    input: "docs/tables/meta_snps_full_eur.cojo.txt"
    output: "results/finemapping/genomewide_loci_{cohorts}_{ancestries}_hg19_v{version}.rp.credible"
    log: "logs/finemapping/genomewide_loci_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    conda: "../envs/finemapping.yaml"
    shell: "Rscript scripts/finemapping/get_genomewide_loci.R {input} {output}"

rule get_cojo_loci:
    input: "docs/tables/meta_snps_full_eur.cojo.format.txt"
    output: finemapping="results/finemapping/cojo_loci_{cohorts}_{ancestries}_hg19_v{version}.rp",
        credible="results/finemapping/cojo_loci_{cohorts}_{ancestries}_hg19_v{version}.rp.credible"
    log: "logs/finemapping/cojo_loci_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    conda: "../envs/finemapping.yaml"
    shell: "bash scripts/finemapping/get_cojo_loci.bash {input} {output.finemapping} {output.credible}"

rule credible_causal_sets:
    input: finemapping="logs/finemapping/create_{cohorts}_{ancestries}_hg19_v{version}.rp.log",
        loci="results/finemapping/genomewide_loci_{cohorts}_{ancestries}_hg19_v{version}.rp.credible"
    output: "results/finemapping/credible_causal_{cohorts}_{ancestries}_hg19_v{version}.rp.complete"
    params: inprefix="results/finemapping/results_{cohorts}_{ancestries}_hg19_v{version}.rp"
    log: "logs/finemapping/credible_causal_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    conda: "../envs/finemapping.yaml" 
    shell: "bash scripts/finemapping/wrap_credible_causal_sets.bash {input.loci} {params.inprefix}"

rule finemap_loci:
    input: snpvar="logs/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}.rp.log",
        n="results/finemapping/n_{cohorts}_{ancestries}_hg19_v{version}.rp.out",
        susie="logs/finemapping/install_susie.log",
        loci="results/finemapping/cojo_loci_{cohorts}_{ancestries}_hg19_v{version}.rp"
    output: "results/finemapping/locus_results_{cohorts}_{ancestries}_hg19_v{version}.rp"
    log: "logs/finemapping/locus_results_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    params: snpvar="results/finemapping/priors_{cohorts}_{ancestries}_hg19_v{version}.rp_chr",
        outprefix="results/finemapping/locus_results_{cohorts}_{ancestries}_hg19_v{version}.rp"
    conda: "../envs/finemapping.yaml"
    shell: "export PYTHONNOUSERSITE=True && scripts/finemapping/finemap_loci.sh {input.loci} {params.snpvar} {input.n} {params.outprefix}"

rule cojo_credible_causal_sets:
    input: finemapping="logs/finemapping/locus_results_{cohorts}_{ancestries}_hg19_v{version}.rp.log",
        loci="results/finemapping/cojo_loci_{cohorts}_{ancestries}_hg19_v{version}.rp.credible"
    output: "results/finemapping/locus_credible_causal_{cohorts}_{ancestries}_hg19_v{version}.rp.complete"
    params: inprefix="results/finemapping/locus_results_{cohorts}_{ancestries}_hg19_v{version}.rp"
    log: "logs/finemapping/locus_credible_causal_{cohorts}_{ancestries}_hg19_v{version}.rp.log"
    conda: "../envs/finemapping.yaml" 
    shell: "bash scripts/finemapping/wrap_locus_credible_causal_sets.bash {input.loci} {params.inprefix}"

