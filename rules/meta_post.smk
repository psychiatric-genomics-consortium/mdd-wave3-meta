############################################
## Post processing of meta-analysed sumstats
############################################

###
### Post processing
###

# check final effect sizes/directions against reference
rule postimp_ma:
    input: "results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.gz"
    output: "results/meta/ma/pgc_mdd_{analysis}.ma"
    shell: """zcat {input} | awk '{{if(NR == 1) {{print $0}} else {{print $0 | "sort -k1,1n -k3,3n"}}}}' | awk '{{if(NR == 1) {{print "SNP", "A1", "A2", "freq", "beta", "se", "p", "N"}} else {{print $2, $4, $5, $7, log($9), $10, $11, 2*$19}}}}' | awk 'NR == 1 || ($4 != 0 && $8 >= 600000)' > {output}"""
    
rule dentist_install:
    input: HTTP.remote("https://www.dropbox.com/s/1mtskir8qzqsmee/DENTIST.1.1.0.0.gz?dl=0", keep_local=False)
    output: "resources/meta/DENTIST_1.1.0.0"
    shell: "zcat {input} > {output}"
    
# SNPs that are not monomorphic in the reference panel
rule postimp_dentist_nomono:
    output: "results/meta/dentist/snps/{ancestries}.chr{chr}.snplist"
    params: popname=lambda wildcards: wildcards.ancestries.upper()
    shell: "zcat {config[refdir]}/HRC.r1-1.EGA.GRCh37.chr{wildcards.chr}.impute.plink.{params.popname}.frq2.gz | awk 'NR != 1 && $6 != 0 {{print $1}}' > {output} "

rule postimp_dentist:
    input: ma="results/meta/ma/pgc_mdd_{cohort}_{ancestries}_hg19_v{version}.ma", snplist="results/meta/dentist/snps/{ancestries}.chr{chr}.snplist"
    params: popname=lambda wildcards: wildcards.ancestries.upper(), prefix="results/meta/dentist/{cohort}_{ancestries}_v{version}/pgc_mdd_{cohort}_{ancestries}_hg19_v{version}.chr{chr}"
    output: "results/meta/dentist/{cohort}_{ancestries}_v{version}/pgc_mdd_{cohort}_{ancestries}_hg19_v{version}.chr{chr}.DENTIST.txt"
    shell: "resources/meta/DENTIST_1.1.0.0 --gwas-summary {input.ma} --bfile {config[refdir]}/pop_{params.popname}/HRC.r1-1.EGA.GRCh37.chr{wildcards.chr}.impute.plink.{params.popname} --chrID {wildcards.chr} --extract {input.snplist} --maf 0.001 --delta-MAF 0.15 --out {params.prefix} --thread-num 16"
    
rule postimp_dentist_chr:
    input: expand("results/meta/dentist/{cohort}_{ancestries}_v{version}/pgc_mdd_{cohort}_{ancestries}_hg19_v{version}.chr{chr}.DENTIST.txt", cohort='full', ancestries='eur', version=analysis_version, chr=range(2, 22))
    output: "results/meta/dentist/pgc_mdd_{cohort}_{ancestries}_hg19_v{version}.DENTIST.txt.gz"
    shell: "cat {input} | gzip -c > {output}"

# check Ricopili output for complete rows and duplicates
rule postimp_rp:
    input: autosome="results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.gz", xsome="results/meta/X/report_pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.gz"
    log: "logs/meta/distribution/{analysis}.rp.log"
    conda: "../envs/meta.yaml"
    output: "results/meta/distribution/pgc_mdd_{analysis}/daner_pgc_mdd_{analysis}.rp.gz"
    script: "../scripts/meta/rp.R"
    
# inputs for postimp_rp
rule postimp_rp_all:
    input: expand("results/meta/distribution/pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.rp.gz", cohorts=cohorts_analyst, ancestries=['eur'], version=analysis_version)