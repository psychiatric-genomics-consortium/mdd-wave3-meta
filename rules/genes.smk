# Gene prioritisation

# install PsyOPS https://github.com/Wainberg/PsyOPS
rule genes_install_psyops:
    output: directory("resources/genes/PsyOPS")
    shell: "git clone https://github.com/Wainberg/PsyOPS.git resources/genes/PsyOPS"

# get top SNPs
rule genes_top_snps:
    input: clump="docs/tables/meta_snps_full_eur.clump.txt", cojo="docs/tables/meta_snps_full_eur.cojo.txt"
    output: clump="results/genes/psyops/meta_snps_full_eur.clump.tsv", cojo="results/genes/psyops/meta_snps_full_eur.cojo.tsv"
    shell: """
    cat {input.clump} | awk -v OFS='\\t' '{{if(NR == 1) {{print "chrom", "bp_hg19", "rs"}} else {{print "chr"$2, $3, $1}}}}' > {output.clump};
    cat {input.cojo} | awk -v OFS='\\t' '{{if(NR == 1) {{print "chrom", "bp_hg19", "rs"}} else if($2 == 1) {{print "chr"$3, $5, $4}}}}' > {output.cojo};
    """

rule genes_psyops:
    input: hits="results/genes/psyops/meta_snps_{cohorts}_{ancestries}.{post}.tsv", psyops=ancient(rules.genes_install_psyops.output)
    output: genes="docs/tables/psyops_{cohorts}_{ancestries}.{post}.txt", log="docs/objects/psyops_{cohorts}_{ancestries}.{post}.log"
    conda: "../envs/genes.yaml"
    shell: """
    cp {input.hits} {input.psyops}/;
    cd {input.psyops}; python PsyOPS.py --GWAS-hit-file $(basename {input.hits}) --output-file $(basename {output.genes}) | tee > $(basename {output.log});
    cd -;
    mv {input.psyops}/$(basename {output.genes}) {output.genes};
    mv {input.psyops}/$(basename {output.log}) {output.log} 
    """
rule genes_psyops_full:
    input: "docs/tables/psyops_full_eur.clump.txt", "docs/tables/psyops_full_eur.cojo.txt"