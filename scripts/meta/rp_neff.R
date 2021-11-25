# QC Ricopili output for Neff
library(dplyr)
library(readr)
library(stringr)

# log
log_path <- snakemake@log[[1]]

# read daner file
daner_gz <- snakemake@input[[1]]
daner <- read_tsv(daner_gz, col_types=cols("SNP"=col_character()))

# QC params
qc_neff <- snakemake@params$qc_neff
qc_info <- snakemake@params$qc_info

# get max half effective sample size for autosomes and X chromosome
daner %>%
mutate(some=case_when(CHR %in% 1:22 ~ 'auto',
                      CHR == 23 ~ 'allo',
                      TRUE ~ NA_character_)) %>%
group_by(some) %>%
summarize(max_Neff_half=max(Neff_half))

daner_neff <-
daner %>% 
mutate(some=case_when(CHR %in% 1:22 ~ 'auto',
       CHR == 23 ~ 'allo',
       TRUE ~ NA_character_)) %>%
group_by(some) %>%
filter(Neff_half >= max(Neff_half)*qc_neff) %>%
ungroup() %>%
select(-some) %>%
mutate(Neff=2*Neff_half)

daner_out_gz <- snakemake@output[[1]]
write_tsv(daner_neff, daner_out_gz)

log_info <- paste("Ricopili Neff:", "remove", nrow(daner)-nrow(daner_neff), "rows with Neff/2 <", 100*qc_neff, "%, keep", nrow(daner_neff), "rows\n")

cat(log_info)

cat(log_info, file=log_path)