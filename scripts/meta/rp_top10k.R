# Find top10k SNPs

library(readr)
library(dplyr)

# top re-clumped SNPs
auto <- read_tsv(snakemake@input$auto)
allo <- read_table(snakemake@input$allo)

# Neff-QC'd SNPs
daner <- read_tsv(snakemake@input$daner)

# COJO and clumped results
cojo <- read_tsv(snakemake@input$cojo)
clump <- read_tsv(snakemake@input$clump)

# merge top SNPs together and arrange by p-value
auto_allo_topNk <-
bind_rows(select(auto, SNP, P),
          select(allo, SNP, P)) %>%
arrange(P)

# clumped or cojo
clump_or_cojo <- unique(c(cojo$SNP, clump$SNP))

# other topN SNPs that were not in the cojo/clumped list
other_top10ke <-
auto_allo_topNk %>%
filter(!SNP %in% clump_or_cojo) %>%
slice(1:(10000-length(clump_or_cojo)))

daner_10k <-
daner %>% filter(SNP %in% c(clump_or_cojo, other_top10ke$SNP)) %>%
mutate(Neff=Neff_half*2)

write_tsv(daner_10k, snakemake@output[[1]])