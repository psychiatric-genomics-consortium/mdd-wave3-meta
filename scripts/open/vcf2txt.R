library(readr)
library(dplyr)
library(gwasvcf)

# load the VCF
vcf <- readVcf(snakemake@input$vcf)

# load the list of HM3 SNPs
hm3 <- read_table2(snakemake@input$hm3)

# gwas info
gwasinfo <- read_tsv(snakemake@input$gwasinfo)

# the dataset we are processing
dataset <- snakemake@wildcards$dataset

# get sample size from GWAS info
sample_size <- gwasinfo %>% filter(id==dataset) %>% pull(sample_size)

# subset to HM3 SNPs
vcf_hm3 <- query_gwas(vcf, rsid=hm3$SNP)

sumstats <- 
vcf_to_granges(vcf_hm3) %>% as_tibble() %>%
transmute(CHR=seqnames, SNP=ID, POS=start, A1=ALT, A2=REF,
          BETA=ES, SE=SE, P=10^(-LP),
          N=coalesce(SS, sample_size))

write_tsv(sumstats, snakemake@output[[1]])