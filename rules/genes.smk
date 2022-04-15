####################
## Gene Priorisation
####################

## GeneMatrix.tsv.gz downloaded from https://figshare.com/articles/dataset/geneMatrix/13335548 on 31st January 2022
rule genes_genematrix:
    output: "resources/fastBAT/geneMatrix.tsv.gz"
    shell: "curl -L https://figshare.com/ndownloader/files/28235607 > {output}"


## Remove double quotes, select only protein_coding genes (N = 19878) and those with known positions on hg19
rule genes_genematrix_protein:
    input: "resources/fastBAT/geneMatrix.tsv.gz"
    output: "resources/fastBAT/geneMatrix.protein.txt"
    shell: """
    zcat {input} | sed 's/\"//g' | awk '($5 == "protein_coding" && $12 ~ /^chr/)' | awk -v range=3 '{{print substr($12,range+1),$13,$14,$3}}' > {output}
    """

## Remove double quotes, select all genes (N = 58844) with known positions on hg19
rule genes_genematrix_all:
    input: "resources/fastBAT/geneMatrix.tsv.gz"
    output: "resources/fastBAT/geneMatrix.all.txt"
    shell: """
    zcat {input} | sed 's/\"//g' | awk '($12 ~ /^chr/)' | awk -v range=3 '{{print substr($12,range+1),$13,$14,$3}}' > {output}
    """

## Extract header row and apply MAF >= 0.01 and infoscore > 0.8 to summary stats
rule genes_fastbat_assoc:
    input: "results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.neff.gz"
    output: "results/fastBAT/pgc_mdd_{cohorts}_{ancestries}_v{version}.assoc.txt"
    shell: "zcat {input} | awk '{{if(NR == 1) {{print $2, $11}} else if($7 >= 0.01 && $7 <= 0.99 && $8 >= 0.8) {{print $2, $11}}}}' > {output}"


## Run fastBAT within GCTA
rule genes_fastbat:
    input: assoc="results/fastBAT/pgc_mdd_{cohorts}_{ancestries}_v{version}.assoc.txt", genelist="resources/fastBAT/geneMatrix.{genes}.txt", bed="resources/1kg/all_phase3.{ancestries}.bed"
    output: "results/fastBAT/pgc_mdd_{cohorts}_{ancestries}_v{version}.{genes}.fbat"
    params: prefix="results/fastBAT/pgc_mdd_{cohorts}_{ancestries}_v{version}.{genes}", bfile="resources/1kg/all_phase3.{ancestries}"
    conda: "../envs/genes.yaml"
    shell: """gcta64 \
    --bfile {params.bfile} \
    --fastBAT {input.assoc} \
    --fastBAT-gene-list {input.genelist} \
    --out {params.prefix} \
    --threads 4
    """

# Run fastbat geneMatrix protein coding gene list (N = 19878) and all geneMatrix gene list (N = 58844)
rule genes_fastbat_analyse:
    input: expand("results/fastBAT/pgc_mdd_full_{ancestries}_v{version}.{genes}.fbat", ancestries=['eur'], version=analysis_version, genes=['protein', 'all'])
