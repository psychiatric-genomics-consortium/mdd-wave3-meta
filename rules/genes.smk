
######################
## Gene Prioritisation
######################

##
## FastBAT
##

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

## Run fastBAT within GCTA using qc'd summary stats and geneMatrix protein coding gene list (N = 19878)
    
## Run fastBAT within GCTA using qc'd summary stats and all geneMatrix gene list (N = 58844)
rule genes_fastbat_analyse:
    input: expand("results/fastBAT/pgc_mdd_full_{ancestries}_v{version}.{genes}.fbat", ancestries=['eur'], version=analysis_version, genes=['protein', 'all'])

##
## PysOPS
##

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

##
## S-LDSC
##

# Make annotation files for partitioned heritability

## Drug Targetor

# Make S-LDSC annotations for drug/genes in the DrugTargetor database
rule genes_sldsc_drugtargetor_annot:
    input: drugtargetor= "resources/drug_enrichment/wholedatabase_for_targetor", geneMatrix="resources/fastBAT/geneMatrix.tsv.gz", phase3="resources/ldsc/1000G_EUR_Phase3_plink"
    params: windowsize=100000, bim="resources/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}.bim"
    output: "resources/drug_enrichment/sldsc/targetor_whole.{chr}.annot.gz"
    conda: "../envs/meta.yaml"
    script: "../scripts/genes/drugtargetor_annot.R"
 
# Compute LD scores with DrugTargettor annotion
rule genes_sldsc_drugtargetor_annot_l2:
    input: annot="resources/drug_enrichment/sldsc/targetor_whole.{chr}.annot.gz",  phase3="resources/ldsc/1000G_EUR_Phase3_plink", ldsc="resources/ldsc/ldsc",  hapmap="resources/ldsc/hapmap3_snps"
    params: bfile="resources/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}", hm="resources/ldsc/hapmap3_snps/hm.{chr}.snp", prefix="resources/drug_enrichment/sldsc/targetor_whole.{chr}"
    output: "resources/drug_enrichment/sldsc/targetor_whole.{chr}.l2.ldscore.gz"
    conda: "../envs/ldsc.yaml"
    shell: """
    python {input.ldsc}/ldsc.py \
    --l2 \
    --bfile {params.bfile} \
    --ld-wind-cm 1 \
    --annot {input.annot} \
    --out {params.prefix} \
    --print-snps {params.hm}
    """
    
rule genes_sldsc_drugtargetor_annot_l2_chr:
    input: expand("resources/drug_enrichment/sldsc/targetor_whole.{chr}.l2.ldscore.gz", chr=range(1, 23))

# GO pathways from FUMA. Test pathways identified by MAGMA with S-LDSC
# list of pathways with Ensemble IDs in results/fuma/go_genesets
rule genes_sldsc_fuma_go_annot:
	input: geneset="results/fuma/go_genesets/{geneset}.GeneSet", phase3="resources/ldsc/1000G_EUR_Phase3_plink", ldsc="resources/ldsc/ldsc",  hapmap="resources/ldsc/hapmap3_snps", gene_coord="resources/ldsc/ENSG_coord.txt"
	params: bim="resources/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}.bim"
	output: "results/fuma/sldsc/{geneset}/{geneset}.{chr}.annot.gz"
	conda: "../envs/ldsc.yaml"
	shell: """
	python {input.ldsc}/make_annot.py \
	--gene-set-file {input.geneset} \
	--gene-coord-file {input.gene_coord} \
	--windowsize 100000 \
	--bimfile {params.bim} \
	--annot-file {output}
	"""
	
# Estimate ld scores
rule genes_sldsc_fuma_go_l2:
	input: annot= "results/fuma/sldsc/{geneset}/{geneset}.{chr}.annot.gz",phase3="resources/ldsc/1000G_EUR_Phase3_plink",hapmap="resources/ldsc/hapmap3_snps", ldsc="resources/ldsc/ldsc"
	params: bfile="resources/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}", hm="resources/ldsc/hapmap3_snps/hm.{chr}.snp", prefix="results/fuma/sldsc/{geneset}/{geneset}.{chr}"
	output: l2="results/fuma/sldsc/{geneset}/{geneset}.{chr}.l2.ldscore.gz", M="results/fuma/sldsc/{geneset}/{geneset}.{chr}.l2.M_5_50"
	conda: "../envs/ldsc.yaml"
	shell: """
	python {input.ldsc}/ldsc.py \
	--l2 \
	--bfile {params.bfile} \
	--ld-wind-cm 1 \
	--annot {input.annot} \
	--thin-annot \
	--out {params.prefix} \
	--print-snps {params.hm}
	"""
	
genes_ldsc_fuma_genesets, = glob_wildcards("results/fuma/go_genesets/{geneset}.GeneSet")

# list pathway/chromsome pairs to generate all files
rule genes_ldsc_fuma_go_l2_all:
	input: expand("results/fuma/sldsc/{geneset}/{geneset}.{chr}.l2.ldscore.gz", geneset=genes_ldsc_fuma_genesets, chr=range(1, 23))