##
## Reference files LDSC and installation
##

# fetch LDSC reference files
rule ldsc_fetch_hm3_bz:
	input: HTTP.remote("https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2")
	output: "resources/ldsc/w_hm3.snplist.bz2"
	shell: "cp {input} {output}"
	
rule ldsc_unzip_hm3:
	input: "resources/ldsc/w_hm3.snplist.bz2"
	output: "resources/ldsc/w_hm3.snplist"
	shell: "bunzip2 {input}"

rule ldsc_fetch_eur_w_ld_chr_bz:
	input: HTTP.remote("https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2")
	output: "resources/ldsc/eur_w_ld_chr.tar.bz2"
	shell: "cp {input} {output}"

rule ldsc_unzip_eur_w_ld_chr:
	input: "resources/ldsc/eur_w_ld_chr.tar.bz2"
	output: directory("resources/ldsc/eur_w_ld_chr")
	shell: "tar -jxvf {input} -C $(dirname {output})"

# LDScore files from PanUKBB: https://pan.ukbb.broadinstitute.org/downloads/index.html for AFR, AMR, CSA, EAS, EUR, and MID
rule ldsc_fetch_ldscore_tar_gz:
    input: HTTP.remote("https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz")
    output: "resources/ldsc/UKBB.ALL.ldscore.tar.gz"
    shell: "cp {input} {output}"

rule ldsc_untar_ldscore:
    input: "resources/ldsc/UKBB.ALL.ldscore.tar.gz"
    output: dir=directory("resources/ldsc/UKBB.ALL.ldscore"), ldscore=expand("resources/ldsc/UKBB.ALL.ldscore/UKBB.{ancestries}.{ext}", ancestries=['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID'], ext=['l2.M', 'l2.M_5_50', 'l2.ldscore.gz', 'rsid.l2.ldscore.gz'])
    shell: "tar -jxvf {input} -C $(dirname {output.dir})"

# install LDSC from github
rule ldsc_install:
	input: "resources/ldsc/w_hm3.snplist", rules.ldsc_unzip_eur_w_ld_chr.output
	output: directory("resources/ldsc/ldsc")
	shell: "git clone https://github.com/bulik/ldsc.git {output}"

##
## LDSC munging and heritability
##
	
# Run munge with N taken from Neff (doi:10.1101/2021.09.22.21263909), FRQ from FRQ_U column
rule ldsc_munge:
	input: sumstats="results/distribution/daner_{cohort}.neff.gz", hm3="resources/ldsc/w_hm3.snplist", ldsc=rules.ldsc_install.output
	params:
		prefix="results/ldsc/munged/{cohort}"
	conda: "../envs/ldsc.yaml"
	output: "results/ldsc/munged/{cohort}.sumstats.gz"
	shell: "resources/ldsc/ldsc/munge_sumstats.py --sumstats {input.sumstats} --snp SNP --N-col Neff --signed-sumstats OR,1 --frq $(gunzip -c {input.sumstats} | head -n 1 | awk '{{print $7}}') --out {params.prefix} --merge-alleles {input.hm3} --chunksize 500000"	
	
# estimate LDscore. Use rsid version of ldscore file and manually cat out M since the file has
# a different prefix
rule ldsc_h2:
	input: sumstats="results/ldsc/munged/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.sumstats.gz", ld=lambda wildcards: expand('resources/ldsc/UKBB.ALL.ldscore/UKBB.{anc}.rsid.l2.ldscore.gz', anc=wildcards.ancestries.upper()), l2_M=lambda wildcards: expand('resources/ldsc/UKBB.ALL.ldscore/UKBB.{anc}.l2.M_5_50', anc=wildcards.ancestries.upper())
	params:
		prefix="results/ldsc/h2/pgc_mdd_{cohorts}_{ancestries}_v{version}", ld=lambda wildcards: expand('resources/ldsc/UKBB.ALL.ldscore/UKBB.{anc}.rsid', anc=wildcards.ancestries.upper())
	conda: "../envs/ldsc.yaml"
	output: "results/ldsc/h2/pgc_mdd_{cohorts}_{ancestries}_v{version}.log"
	shell: "resources/ldsc/ldsc/ldsc.py --h2 {input.sumstats} --ref-ld {params.ld} --w-ld {params.ld} --M $(cat {input.l2_M}) --out {params.prefix}"


##
## Stratified LDSC resources
##

# 1000G Phase 3 plink files for LDSC
rule ldsc_1kg3:
    input: HTTP.remote("https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_plinkfiles.tgz")
    output: directory("resources/ldsc/1000G_EUR_Phase3_plink")
    shell: "tar xzf {input} -C $(dirname {output})"
    
# Gene coordinates file
rule ldsc_gene_coord:
    input: HTTP.remote("https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/make_annot_sample_files/ENSG_coord.txt")
    output: "resources/ldsc/ENSG_coord.txt"
    shell: "mv {input} {output}"
    
# hapmap3 snps
rule ldsc_hapmap3:
    input: HTTP.remote("https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/hapmap3_snps.tgz")
    output: directory("resources/ldsc/hapmap3_snps")
    shell: "tar xzf {input} -C $(dirname {output})"
    
# regression weights
rule ldsc_weights:
    input: HTTP.remote("https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/weights_hm3_no_hla.tgz")
    output: directory("resources/ldsc/weights_hm3_no_hla")
    shell: "tar xzf {input} -C $(dirname {output})"
    
# allele frequencies
rule ldsc_mac5eur:
    input: HTTP.remote("https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_frq.tgz")
    output: directory("resources/ldsc/1000G_Phase3_frq")
    shell: "tar xzf {input} -C $(dirname {output})"
    
# baseline annotations
rule ldsc_baseline:
    input: HTTP.remote("https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baseline_v1.2_ldscores.tgz")
    output: directory("resources/ldsc/baseline_v1.2")
    shell: "tar xzf {input} -C $(dirname {output})"

# baseline SNPs
rule ldsc_baseline_snps:
	input: "resources/ldsc/baseline_v1.2"
	output: "resources/ldsc/baseline_v1.2_snps/baseline.{chr}.snp"
	shell: "gunzip -c resources/ldsc/baseline_v1.2/baseline.{wildcards.chr}.l2.ldscore.gz | awk 'NR > 1 {{print $2}}' > {output}"