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

# get max half effective sample size
max_neff_half <- max(daner$Neff_half)

daner_neff <-
daner %>% 
filter(Neff_half >= max_neff_half*qc_neff) %>%
mutate(Neff=2*Neff_half)

daner_out_gz <- snakemake@output[[1]]
write_tsv(daner_neff, daner_out_gz)

log_info <- paste("Ricopili Neff:", "remove", nrow(daner)-nrow(daner_neff), "rows with Neff/2 <", max_neff_half*qc_neff, ", keep", nrow(daner_neff), "rows\n")

cat(log_info)

cat(log_info, file=log_path)