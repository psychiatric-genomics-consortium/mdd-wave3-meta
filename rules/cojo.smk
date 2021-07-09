# Conditional analysis

# Use UK Biobank as reference
# add entry `ukb: path/to/ukbdata` directory path to the config.yaml file.
# Directory should contain
# - ukb_imp_chrN_v3.bgen: bgen files
# - ukb_imp_chrN_v3.bgen.bgi: bgi files
# - ukbAAA_imp_chrN_v3_sSSSS.sample: sample files
# - ukbRRRR.enc_ukb: phenotype release file (unencrypted)


# download qctool
rule cojo_qctool_tgz:
    input: HTTP.remote("https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz", keep_local=False)
    output: "resources/cojo/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz"
    shell: "mv {input} {output}"

rule cojo_qctool:
    input: "resources/cojo/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz"
    output: "resources/cojo/qctool_v2.0.8-CentOS\\ Linux7.6.1810-x86_64/qctool"
    shell: "tar -xz --directory=resources/cojo -f {input}"

# QC sumstats to MAF <= 0.01, INFO >= 0.6
rule cojo_daner_qc:
    input: "results/distribution/daner_{analysis}.rp.gz"
    output: "results/cojo/daner_{analysis}.rp.qc.gz"
    shell: "zcat {input} | awk '{{if(NR == 1 || ($7 >= 0.01 && $7 <= 0.99 && $8 >= 0.6)) {{print $0}}}}' | gzip -c > {output}"
    
# sumstats for input into GCTA
rule cojo_ma:
    input: "results/cojo/daner_{analysis}.rp.qc.gz"
    output: "results/cojo/{analysis}.ma"
    shell: """zcat {input} | awk '{{if(NR == 1) {{print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"}} else {{print $2, $4, $5, $7, log($9), $10, $11, 2*$19}}}}' > {output}"""
    
# regions to analyse
# Pull regions clumped from Ricopili
rule cojo_regions:
    input: "results/distribution/daner_{analysis}.gz.p4.clump.areator.sorted.1mhc"
    output: "results/cojo/{analysis}.regions"
    shell: "cat {input} | awk 'NR > 1 && $4 <= 5e-8 {{print $2, $14, $15}}' > {output}"

# symlink to UKB bgen and sample files
rule cojo_ukb_bgen:
    input: bgen=expand("{ukb}/ukb_imp_chr{{chr}}_v3.bgen", ukb=config['ukb'])
    output: "resources/cojo/ukb/ukb_imp_chr{chr}_v3.bgen"
    shell: "ln -sv {input} {output}"

rule cojo_ukb_bgi:
    input: bgi=expand("{ukb}/ukb_imp_chr{{chr}}_v3.bgen.bgi", ukb=config['ukb'])
    output: "resources/cojo/ukb/ukb_imp_chr{chr}_v3.bgen.bgi"
    shell: "ln -sv {input} {output}"
    
rule cojo_ukb_sample:
    params: ukb=config['ukb']
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
    params: ukb=config['ukb'], prefix="resources/cojo/ukb/f22006"
    output: "results/cojo/f22006.txt"
    shell: "{input.ukbconv} {params.ukb}/ukb*.enc_ukb txt -o{params.prefix} -s22006"
    
rule cojo_ukb_eur_ids:
    input: "results/cojo/f22006.txt"
    output: "results/cojo/ukb_eur.ids"
    shell: "cat {input} | awk '{{if(NR > 1 && $2 == 1) {{print $1, $1}}}}' > {output}"
    
# extract list related individuals (ukbA_rel_sP.txt from `gfetch rel`)
rule cojo_ukb_rels:
    params: ukb=config['ukb']
    output: "results/cojo/ukb_rel.dat"
    shell: "cat {params.ukb}/ukb*_rel_s*.dat | awk 'NR > 1 {{print $1}}' | sed 's/ /\\n/' | sort | uniq | awk '{{print $1, $1}}' > {output}"
    
# read in list of regions (CHR START STOP) and 
# merge close by regions
# return region strings "CHR:START-STOP"
def cojo_parse_regions(regions_file):
    # open regions file and split into list of [CHR, START, STOP]
    regions_f = open(regions_file, 'r')
    regions_list = regions_f.read().split('\n')
    regions = [[int(float(i)) for i in r.split()] for r in regions_list if r != '']
    # merge together regions that are on the same chromosome and within 50kb of each
    # other
    merged_regions = []
    kb50 = 50000
    for r in regions:
        add = True
        for i in range(len(merged_regions)):
            m = merged_regions[i]
            if r[0] == m[0]:
                if abs(m[1] - r[2]) <= kb50 or abs(r[1] - m[2]) <= kb50:
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
    conda: "../envs/cojo.yaml"
    shell: "bgenix -g {input.bgen} -incl-range {params.chr0}:{wildcards.start}-{wildcards.stop} > {output}"
    
# Extract SNPs for each region
rule cojo_snplists:
    input: "results/distribution/daner_{analysis}.rp.gz"
    output: "results/cojo/{analysis}/{chr}:{start}-{stop}.snplist"
    shell: """zcat {input} | awk '{{if(NR > 1 && $1 == {wildcards.chr} && {wildcards.start} <= $3 && $3 <= {wildcards.stop}) {{print $2}}}}' > {output}"""
    
    
# convert UKB BGEN to PLINK BED
# keep unrelated European ancestries
# extract SNPs that are in the GWAS
rule cojo_region_bed:
    input: bgen="results/cojo/{analysis}/{chr}:{start}-{stop}.bgen", snplist="results/cojo/{analysis}/{chr}:{start}-{stop}.snplist", sample= "resources/cojo/ukb/ukb_imp_chr{chr}_v3.sample", eur_ids="results/cojo/ukb_eur.ids", rel_ids="results/cojo/ukb_rel.dat"
    conda: "../envs/cojo.yaml"
    params: prefix="results/cojo/{analysis}/{chr}:{start}-{stop}"
    output: temp("results/cojo/{analysis}/{chr}:{start}-{stop}.bed"), temp("results/cojo/{analysis}/{chr}:{start}-{stop}.fam"), temp("results/cojo/{analysis}/{chr}:{start}-{stop}.bim")
    shell: "plink2 --make-bed --bgen {input.bgen} 'ref-first' --sample {input.sample} --double-id --keep {input.eur_ids} --remove {input.rel_ids} --extract {input.snplist} --out {params.prefix} --memory 4000 --threads 1"
    
# COJO analysis
rule cojo_slct:
    input: bed="results/cojo/{analysis}/{chr}:{start}-{stop}.bed", fam="results/cojo/{analysis}/{chr}:{start}-{stop}.fam", bim="results/cojo/{analysis}/{chr}:{start}-{stop}.bim", ma="results/cojo/{analysis}.ma"
    conda: "../envs/cojo.yaml"
    params: prefix="results/cojo/{analysis}/{chr}:{start}-{stop}"
    output: jma="results/cojo/{analysis}/{chr}:{start}-{stop}.jma.cojo", jma_ldr="results/cojo/{analysis}/{chr}:{start}-{stop}.ldr.cojo"
    shell: "gcta64 --bfile {params.prefix} --cojo-file {input.ma} --cojo-slct --out {params.prefix}"
    

rule cojo_regions_analyse:
    input: lambda wildcards: expand("results/cojo/{analysis}/{region}.jma.cojo", analysis=wildcards.analysis, region=cojo_parse_regions("results/cojo/" + wildcards.analysis + '.regions'))
    output: "results/cojo/{analysis}.cojo"
    shell: "touch {output}"

rule cojo_analyse:
    input: expand("results/cojo/pgc_mdd_full_eur_hg19_v{version}.cojo", version=analysis_version)