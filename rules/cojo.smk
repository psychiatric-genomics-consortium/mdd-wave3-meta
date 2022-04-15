# Conditional analysis

# Use UK Biobank as reference
# add entry `ukb: path/to/ukbdata` directory path to the config.yaml file.
# Directory should contain
# - ukb_imp_chrN_v3.bgen: bgen files
# - ukb_imp_chrN_v3.bgen.bgi: bgi files
# - ukbAAA_imp_chrN_v3_sSSSS.sample: sample files
# - ukbRRRR.enc_ukb: phenotype release file (unencrypted)

ukb_path=config['ukb'] if "ukb" in config else '~'


# download qctool
rule cojo_qctool_tgz:
    input: HTTP.remote("https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz", keep_local=False)
    output: "resources/cojo/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz"
    shell: "mv {input} {output}"

rule cojo_qctool:
    input: "resources/cojo/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz"
    output: "resources/cojo/qctool_v2.0.8-CentOS\\ Linux7.6.1810-x86_64/qctool"
    shell: "tar -xz --directory=resources/cojo -f {input}"

# QC sumstats to MAF <= 0.01, INFO >= 0.6 based on Neff > 80% max version
rule cojo_daner_qc:
    input: "results/distribution/daner_{analysis}.neff.gz"
    output: "results/cojo/daner_{analysis}.qc.gz"
    params: qc_neff=0.8
    shell: "zcat {input} | awk '{{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99 && $8 >= 0.6)) {{print $0}}}}' | gzip -c > {output}"
    
# sumstats for input into GCTA
rule cojo_ma:
    input: "results/cojo/daner_{analysis}.qc.gz"
    output: "results/cojo/{analysis}.ma"
    conda: "../envs/meta.yaml"
    script: "../scripts/meta/cojo_ma.R"
    
# regions to analyse
# Pull regions clumped from Ricopili
rule cojo_regions:
    input: "results/distribution/daner_{analysis}.gz.p4.clump.areator.sorted.1mhc"
    output: "results/cojo/{analysis}.regions"
    shell: "cat {input} | awk 'NR > 1 && $4 <= 5e-8 {{print $2, $14, $15}}' > {output}"

# symlink to UKB bgen and sample files
rule cojo_ukb_bgen:
    input: bgen=expand("{ukb}/ukb_imp_chr{{chr}}_v3.bgen", ukb=ukb_path)
    output: "resources/cojo/ukb/ukb_imp_chr{chr}_v3.bgen"
    shell: "ln -sv {input} {output}"

rule cojo_ukb_bgi:
    input: bgi=expand("{ukb}/ukb_imp_chr{{chr}}_v3.bgen.bgi", ukb=ukb_path)
    output: "resources/cojo/ukb/ukb_imp_chr{chr}_v3.bgen.bgi"
    shell: "ln -sv {input} {output}"
    
rule cojo_ukb_sample:
    params: ukb=ukb_path
    output: "resources/cojo/ukb/ukb_imp_chr{chr}_v3.sample"
    shell: "ln -sv {params.ukb}/ukb*_imp_chr{wildcards.chr}_v3_s*.sample {output}"

rule cojo_ukb_conv:
    input: HTTP.remote("https://biobank.ndph.ox.ac.uk/ukb/util/ukbconv", keep_local=False)
    output: "resources/cojo/ukb/ukbconv"
    shell: "mv {input} {output}; chmod a+x {output}"

# extract variable identifying European ancestries
# Assumes there is one file in the config['ukb'] directory with the extension *.enc_ukb
rule cojo_ukb_22006:
    input: ukbconv="resources/cojo/ukb/ukbconv"
    params: ukb=ukb_path, prefix="resources/cojo/ukb/f22006"
    output: "results/cojo/f22006.txt"
    shell: "{input.ukbconv} {params.ukb}/ukb*.enc_ukb txt -o{params.prefix} -s22006"
    
rule cojo_ukb_eur_ids:
    input: "results/cojo/f22006.txt"
    output: "results/cojo/ukb_eur.ids"
    shell: "cat {input} | awk '{{if(NR > 1 && $2 == 1) {{print $1, $1}}}}' > {output}"
    
# extract list related individuals (ukbA_rel_sP.txt from `gfetch rel`)
rule cojo_ukb_rels:
    params: ukb=ukb_path
    output: "results/cojo/ukb_rel.dat"
    shell: "cat {params.ukb}/ukb*_rel_s*.dat | awk 'NR > 1 {{print $1}}' | sed 's/ /\\n/' | sort | uniq | awk '{{print $1, $1}}' > {output}"
    
# read in list of regions (CHR START STOP) and 
# merge close by regions
# return region strings "CHR:START-STOP"
def cojo_parse_regions(regions_file, kb=50):
    from operator import itemgetter
    # open regions file and split into list of [CHR, START, STOP]
    regions_f = open(regions_file, 'r')
    regions_list = regions_f.read().split('\n')
    regions = [[int(float(i)) for i in r.split()] for r in regions_list if r != '']
    # merge together regions that are on the same chromosome and within M kb of each
    # other. Sort regions first
    merged_regions = []
    bp = kb*1000
    for r in sorted(regions, key=itemgetter(0, 1)):
        add = True
        for i in range(len(merged_regions)):
            m = merged_regions[i]
            if r[0] == m[0]:
                # extend regions
                r_start=r[1] - bp/2
                r_end=r[2] + bp/2
                m_start=m[1] - bp/2
                m_end=m[2] + bp/2
                # test for partial or full overlap
                if(r_start <= m_end and m_start <= r_end):
                    # merge regions
                    merged_regions[i] = [r[0], min(r[1], m[1]), max(r[2], m[2])]
                    add = False
                    break
            
        if add:
            merged_regions = merged_regions + [r]
            
    region_strings = [str(r[0]) + ':' + str(r[1]) + '-' + str(r[2]) for r in merged_regions]
    return region_strings

# Extract region from the BGEN file
rule cojo_region_bgen:
    input: bgen="resources/cojo/ukb/ukb_imp_chr{chr}_v3.bgen", bgi="resources/cojo/ukb/ukb_imp_chr{chr}_v3.bgen.bgi", regions="results/cojo/{analysis}.regions"
    output: temp("results/cojo/{analysis}/{chr}:{start}-{stop}.bgen")
    params: chr0=lambda wildcards: wildcards.chr.zfill(2)
    group: "cojo"
    conda: "../envs/cojo.yaml"
    shell: "bgenix -g {input.bgen} -incl-range {params.chr0}:{wildcards.start}-{wildcards.stop} > {output}"
    
# Extract SNPs for each region
# List of CPIDs and SNP names for renaming SNPs
rule cojo_varids:
    input: "results/cojo/daner_{analysis}.qc.gz"
    output: "results/cojo/{analysis}/{chr}:{start}-{stop}.varids"
    group: "cojo"
    shell: """zcat {input} | awk '{{if(NR > 1 && $1 == {wildcards.chr} && {wildcards.start} <= $3 && $3 <= {wildcards.stop}) {{print $1, $2, 0, $3, $4, $5}}}}' > {output}"""
 
# list of SNPs to keep 
rule cojo_snplist:
    input: "results/cojo/{analysis}/{chr}:{start}-{stop}.varids"
    output: "results/cojo/{analysis}/{chr}:{start}-{stop}.snplist"
    group: "cojo"
    shell: "cat {input} | awk '{{print $2}}' > {output}"
    
    
# convert UKB BGEN to PLINK BED
# keep unrelated European ancestries
# extract SNPs in the specified region
# rename variant IDs based on sumstats
rule cojo_region_bed:
    input: bgen="results/cojo/{analysis}/{chr}:{start}-{stop}.bgen", varids="results/cojo/{analysis}/{chr}:{start}-{stop}.varids", snplist="results/cojo/{analysis}/{chr}:{start}-{stop}.snplist", sample= "resources/cojo/ukb/ukb_imp_chr{chr}_v3.sample", eur_ids="results/cojo/ukb_eur.ids", rel_ids="results/cojo/ukb_rel.dat"
    conda: "../envs/cojo.yaml"
    params: prefix="results/cojo/{analysis}/{chr}:{start}-{stop}"
    output: temp("results/cojo/{analysis}/{chr}:{start}-{stop}.bed"), temp("results/cojo/{analysis}/{chr}:{start}-{stop}.fam"), "results/cojo/{analysis}/{chr}:{start}-{stop}.bim"
    group: "cojo"
    shell: """plink2 --make-bed --bgen {input.bgen} 'ref-first' \
    --sample {input.sample} --double-id \
    --keep {input.eur_ids} --remove {input.rel_ids} \
    --recover-var-ids {input.varids} 'partial' \
    --extract {input.snplist} \
    --out {params.prefix} --memory 4000 --threads 1"""

# COJO analysis
# parse CHR, start, and stop from input filename
rule cojo_slct:
    input: bed="results/cojo/{analysis}/{chr}:{start}-{stop}.bed", fam="results/cojo/{analysis}/{chr}:{start}-{stop}.fam", bim="results/cojo/{analysis}/{chr}:{start}-{stop}.bim", ma="results/cojo/{analysis}.ma"
    conda: "../envs/cojo.yaml"
    params: prefix="results/cojo/{analysis}/{chr}:{start}-{stop}"
    output: jma="results/cojo/{analysis}/{chr}:{start}-{stop}.jma.cojo", jma_ldr="results/cojo/{analysis}/{chr}:{start}-{stop}.ldr.cojo"
    group: "cojo"
    shell: "gcta64 --bfile {params.prefix} --cojo-file {input.ma} --cojo-slct --out {params.prefix}; if grep -e 'No SNPs have been selected' {params.prefix}.log; then touch {output}; fi"
    
# parse regions from the region list file to determine inputs
# cojo_parse_regions() returns a list ["CHR:START-STOP", "CHR:START-STOP", ...]
rule cojo_regions_analyse:
    input: cojo=lambda wildcards: expand("results/cojo/{analysis}/{region}.jma.cojo", analysis=wildcards.analysis, region=cojo_parse_regions("results/cojo/" + wildcards.analysis + '.regions', meta_qc_params['cojo_kb'])), bim=lambda wildcards: expand("results/cojo/{analysis}/{region}.bim", analysis=wildcards.analysis, region=cojo_parse_regions("results/cojo/" + wildcards.analysis + '.regions', meta_qc_params['cojo_kb'])), daner="results/cojo/daner_{analysis}.qc.gz", clump="results/distribution/daner_{analysis}.gz.p4.clump.areator.sorted.1mhc.txt"
    conda: "../envs/meta.yaml"
    log: "logs/cojo/{analysis}.log"
    output: cojo="results/cojo/{analysis}.cojo", singletons="results/cojo/{analysis}.singleton.cojo"
    script: "../scripts/meta/cojo.R"
    
rule cojo_table_eur:
    input: cojo=expand("results/cojo/pgc_mdd_{{cohorts}}_eur_hg19_v{version}.cojo", version=analysis_version), singletons=expand("results/cojo/pgc_mdd_{{cohorts}}_eur_hg19_v{version}.singleton.cojo", version=analysis_version)
    output: cojo="docs/tables/meta_snps_{cohorts}_eur.cojo.txt", singletons="docs/tables/meta_snps_{cohorts}_eur.cojo.singleton.txt"
    shell: "cp {input.cojo} {output.cojo}; cp {input.singletons} {output.singletons}"

# run all COJO analyses
rule cojo_analyse:
    input: "docs/tables/meta_snps_full_eur.cojo.txt"
    
# copy COJO log into the repository
rule cojo_log:
    input: expand("logs/cojo/pgc_mdd_full_eur_hg19_v{version}.log", version=analysis_version)
    output: "docs/objects/meta_snps_full_eur.cojo.log"
    shell: "cp {input} {output}"
    
# download top SNPs from Wray 2018
rule cojo_wray:
    conda: "../envs/meta.yaml"
    output: "docs/tables/previous/wray2018_table_2.txt"
    shell: """
    Rscript --default-packages=readr,xml2,rvest -e "write_tsv(html_table(read_html('https://www.nature.com/articles/s41588-018-0090-3/tables/2'))[[1]], '{output}')"
    """
    
# download top SNPs from Howard 2019
rule cojo_howard:
    input: HTTP.remote("https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-018-0326-7/MediaObjects/41593_2018_326_MOESM3_ESM.xlsx", keep_local=False)
    output: "docs/tables/previous/howard2019_table_s1.xlsx"
    shell: "cp {input} {output}"
    
# download GWAS catalogue hits for unipolar depression
# https://www.ebi.ac.uk/gwas/efotraits/EFO_0003761
# rule cojo_gwas_catalog:
#     input: HTTP.remote("https://www.ebi.ac.uk/gwas/api/search/downloads/full", keep_local=False)
#     output: "docs/tables/previous/gwas_catalog_EFO_0003761.txt"
#     shell: "cat {input} | grep "

# find range of tagging variants for previous sumstats
rule cojo_previous_variants:
    input: levey="docs/tables/previous/levey2021_223snps.txt", giannakopoulou="docs/tables/previous/Giannakopoulou2021_table.txt", catalog="docs/tables/previous/gwas-association-EFO_0003761-withChildTraits.tsv.bz2", bim="resources/1kg/1kg_phase1_all.bim"
    conda: "../envs/meta.yaml"
    output: "results/cojo/previous/previous.snpid"
    script: "../scripts/meta/cojo_previous.R"
    
rule cojo_previous_tagging_variants:
    input: snpids="results/cojo/previous/previous.snpid", bed="resources/1kg/1kg_phase1_all.bed", bim="resources/1kg/1kg_phase1_all.bim", fam="resources/1kg/1kg_phase1_all.fam", dups="resources/1kg/1kg_phase1_all.dups"
    conda: "../envs/meta.yaml"
    params: prefix="results/cojo/previous/previous"
    output: "results/cojo/previous/previous.tags.list"
    shell: "plink --bfile resources/1kg/1kg_phase1_all --maf 0.01 --exclude {input.dups} --show-tags {input.snpids} --list-all --tag-kb 500 --tag-r2 0.8 --out {params.prefix}"
    

# install genpwr R library
rule cojo_install_genpwr:
    output: touch("resources/cojo/install_genpwr.done")
    conda: "../envs/meta.yaml"
    shell: """
    Rscript -e "devtools::install_github('camillemmoore/Power_Genetics', subdir='genpwr')"
    """

# install ggman R library
rule cojo_install_ggman:
    output: touch("resources/cojo/install_ggman.done")
    conda: "../envs/meta.yaml"
    shell: """
    Rscript -e "devtools::install_github('drveera/ggman')"
    """

rule cojo_docs:
    input: cojo="docs/tables/meta_snps_full_eur.cojo.txt", log="docs/objects/meta_snps_full_eur.cojo.log", wray="docs/tables/previous/wray2018_table_2.txt", howard="docs/tables/previous/howard2019_table_s1.xlsx", levey="docs/tables/previous/levey2021_223snps.txt", giannakopoulou="docs/tables/previous/Giannakopoulou2021_table.txt", catalog="docs/tables/previous/gwas-association-EFO_0003761-withChildTraits.tsv.bz2", rp_clump="docs/tables/meta_snps_full_eur.clump.txt", tags="results/cojo/previous/previous.tags.list", daner=expand("results/distribution/daner_pgc_mdd_full_eur_hg19_v{version}.neff.gz", version=analysis_version), rmd="docs/cojo.Rmd", genpwr=ancient(rules.cojo_install_genpwr.output), ggman=ancient(rules.cojo_install_ggman.output)
    params: qc=meta_qc_params
    output: "docs/cojo.md"
    conda: "../envs/meta.yaml"
    script: "../docs/cojo.Rmd"