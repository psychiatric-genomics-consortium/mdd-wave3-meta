##########
# Run PWAS analysis
##########

####
# Format ROSMAP and Banner PWAS data
####

rule format_pwas_data:
  output:
    "resources/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos"
  conda:
    "../envs/twas.yaml"
  params:
    rosmap_fusion= config["rosmap_fusion"],
    banner_fusion= config["banner_fusion"]
  shell:
    "Rscript scripts/pwas/format_pwas_data.R \
      --rosmap {params.rosmap_fusion} \
      --banner {params.banner_fusion}"

####
# Format ROSMAP SMR data
####

rule format_rosmap_smr_data:
  input:
    rules.download_smr.output,
    rules.prep_1kg.output
  output:
    "resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd.epi"
  conda:
    "../envs/twas.yaml"
  params:
    rosmap_smr= config["rosmap_smr"],
  shell:
    "Rscript scripts/pwas/format_rosmap_smr_data.R \
      --rosmap {params.rosmap_smr}"

###
# Run PWAS
###

# Run twas using ROSMAP SNP-weights
rule run_rosmap_pwas:
  resources: mem_mb=20000 
  input:
    sumstats=rules.focus_munge.output,
    neff_txt="results/twas/median_Neff.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    format_psychencode=rules.format_pwas_data.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "results/pwas/PGC_MDD3_pwas_rosmap_chr{chr}"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Neff=$(cat results/twas/median_Neff.txt); Rscript resources/twas/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights resources/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS/train_weights_withN.pos \
    --weights_dir resources/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS \
    --ref_ld_chr resources/twas/1kg/EUR/EUR_phase3.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 1e-3 \
    --GWASN ${{Neff}}"

rule rosmap_pwas_all_chr:
    input: 
      lambda w: expand("results/pwas/PGC_MDD3_pwas_rosmap_chr{chr}", chr=range(1, 23))
    output: 
      touch("results/pwas/PGC_MDD3_pwas_rosmap_all_chr.done")

# Run twas using Banner SNP-weights
rule run_banner_pwas:
  resources: mem_mb=20000 
  input:
    sumstats=rules.focus_munge.output,
    neff_txt="results/twas/median_Neff.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    format_psychencode=rules.format_pwas_data.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "results/pwas/PGC_MDD3_pwas_banner_chr{chr}"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Neff=$(cat results/twas/median_Neff.txt); Rscript resources/twas/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights resources/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos \
    --weights_dir resources/data/banner_twas/Banner.n152.fusion.WEIGHTS \
    --ref_ld_chr resources/twas/1kg/EUR/EUR_phase3.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 0.05 \
    --GWASN ${{Neff}}"

rule banner_pwas_all_chr:
    input: 
      lambda w: expand("results/pwas/PGC_MDD3_pwas_banner_chr{chr}", chr=range(1, 23))
    output: 
      touch("results/pwas/PGC_MDD3_pwas_banner_all_chr.done")

#####
# Run SMR using ROSMAP
#####

rule run_rosmap_smr:
  input:
    rules.prep_1kg.output,
    rules.daner_to_cojo.output,
    rules.format_rosmap_smr_data.output,
    rules.download_smr.output
  output:
    "results/pwas/rosmap_smr/rosmap_smr_res_chr{chr}.smr"
  conda: 
    "../envs/twas.yaml"
  shell:
    "resources/twas/smr/smr_Linux \
      --bfile resources/twas/1kg/EUR/EUR_phase3.MAF_001.chr{wildcards.chr} \
      --gwas-summary results/twas/munged_gwas/daner_pgc_mdd_full_eur_hg19_v3.49.24.11_COJO.txt \
      --beqtl-summary resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd \
      --out results/pwas/rosmap_smr/rosmap_smr_res_chr{wildcards.chr}"

rule run_rosmap_smr_chr:
    input: 
      lambda w: expand("results/pwas/rosmap_smr/rosmap_smr_res_chr{chr}.smr", chr=range(1, 23))
    output: 
      touch("results/pwas/rosmap_smr/rosmap_smr_res_all_chr.done")

# Format rosmap smr results
rule process_rosmap_smr:
  input:
    "results/pwas/rosmap_smr/rosmap_smr_res_all_chr.done"
  output:
    "results/pwas/rosmap_smr/rosmap_smr_res_GW.txt.gz"
  conda: 
    "../envs/twas.yaml"
  shell:
    "Rscript scripts/pwas/process_rosmap_smr.R"

########
# Create report of the results
########

rule create_pwas_report:
  input:
    rules.process_rosmap_smr.output,
    rules.banner_pwas_all_chr.output,
    rules.rosmap_pwas_all_chr.output
  output:
    "docs/pwas_report.html"
  conda:
    "../envs/twas.yaml"
  shell:
    "Rscript -e \"rmarkdown::render(\'scripts/pwas/pwas_report.Rmd\', output_file = \'../../docs/pwas_report.html\')\""

