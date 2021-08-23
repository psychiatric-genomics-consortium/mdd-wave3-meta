library(gwasvcf)
library(readr)

# load the VCF
vcf <- readVcf(snakemake@input$vcf)

# load the list of HM3 SNPs
hm3 <- read_table2(snakemake@input$hm3)

# subset to HM3 SNPs
vcf_hm3 <- query_gwas(vcf, rsid=hm3$SNP)

sumstats <- 
vcf_to_granges(vcf_hm3) %>% dplyr::as_tibble() %>%
dplyr::transmute(CHR=seqnames, SNP=ID, POS=start, A1=ALT, A2=REF, BETA=ES, SE=SE, N=SS, P=10^(-LP))

write_tsv(sumstats, snakemake@output[[1]])