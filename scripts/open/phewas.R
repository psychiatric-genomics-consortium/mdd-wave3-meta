library(ieugwasr)
library(readr)
library(dplyr)

# read in the top SNPs
cojo <- read_tsv(snakemake@input[[1]])

cojo_region_top <- cojo %>%
group_by(region) %>%
filter(P == min(P))

# batches to search. See https://gwas.mrcieu.ac.uk/datasets/
batch <- c('ebi-a', 'finn-a', 'ieu-a', 'ieu-b', 'met-b', 'ubm-a', 'ukb-a', 'ukb-b',  'ukb-d')

# query for phewas
phewas_assoc <- phewas(variants=cojo_region_top$SNP, batch=batch, pval=1e-5)

write_tsv(phewas_assoc, snakemake@output$phewas)

# query GWAS info
phewas_info <- gwasinfo(unique(phewas_assoc$id))

write_tsv(phewas_info, snakemake@output$gwasinfo)