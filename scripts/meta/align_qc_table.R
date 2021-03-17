library(dplyr)
library(readr)
library(gdata)

qc_align <-
bind_rows(lapply(unlist(snakemake@input), read_tsv, col_types=cols(.default=col_character()))) %>%
arrange(ancestries, cohort, release)

write.fwf(as.data.frame(qc_align), file=snakemake@output[[1]], sep="\t")