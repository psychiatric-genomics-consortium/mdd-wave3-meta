# Find top10k SNPs

library(readr)
library(dplyr)

# top re-clumped SNPs
auto <- read_tsv(snakemake@input$auto)
allo <- read_table(snakemake@input$allo)

# Neff-QC'd SNPs
daner <- read_tsv(snakemake@input$daner)

# COJO and clumped results
cojo_list <- lapply(snakemake@input$cojo, read_tsv)
cojo <- bind_rows(lapply(cojo_list, select, SNP))
clump_list <- lapply(snakemake@input$clump, read_tsv)
clump <- bind_rows(lapply(clump_list, select, SNP))

# merge top SNPs together and arrange by p-value
auto_allo_top_n <-
  bind_rows(select(auto, SNP, P),
            select(allo, SNP, P)) |>
  arrange(P)

# clumped or cojo
clump_or_cojo <- bind_rows(cojo, clump) |>
  distinct() |>
  pull(SNP)

# other topN SNPs that were not in the cojo/clumped list
other_top10ke <-
  auto_allo_top_n |>
  filter(!SNP %in% clump_or_cojo) %>%
  slice(1:(10000 - length(clump_or_cojo)))

daner_10k <- daner %>%
  filter(SNP %in% c(clump_or_cojo, other_top10ke$SNP)) |>
  mutate(Neff = Neff_half * 2)

write_tsv(daner_10k, snakemake@output[[1]])