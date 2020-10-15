library(gwasvcf)
library(VariantAnnotation)
library(dplyr)
library(magrittr)

sumstats_vcf_file <- VcfFile(snakemake@input$sumstats)
sumstats_vcf <- readVcf(sumstats_vcf_file, genome=snakemake@wildcards$build)