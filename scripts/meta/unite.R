library(dplyr)
library(readr)

cpid_files <- snakemake@input

cpids <- lapply(cpid_files, read_tsv)

cpids_counts <- bind_rows(cpids) %>%
count(SNP)