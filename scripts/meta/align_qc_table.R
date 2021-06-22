library(dplyr)
library(readr)
library(gdata)

aligned_txt <- snakemake@input$aligned
filtered_txt <- snakemake@input$filtered

qc_align <-
bind_rows(lapply(unlist(aligned_txt), read_tsv, col_types=cols(.default=col_character()))) %>%
arrange(ancestries, cohort, release)

qc_filter <-
bind_rows(lapply(unlist(filtered_txt), read_tsv, col_types=cols(.default=col_character()))) %>%
arrange(ancestries, cohort, release)

qc <- qc_align %>% 
left_join(qc_filter, by=c('ancestries', 'cohort', 'release', 'N_cases', 'N_controls'))

write.fwf(as.data.frame(qc), file=snakemake@output[[1]], sep="\t")