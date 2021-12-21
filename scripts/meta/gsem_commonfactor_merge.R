library(dplyr)
library(readr)

sumstats <- bind_rows(lapply(snakemake@input$sumstats, read_tsv, col_types=cols(CHR=col_integer(), warning=col_character())))

write_tsv(sumstats, snakemake@output[[1]])