# Check Ricopili output

library(dplyr)
library(readr)
library(stringr)

# log
log_path <- snakemake@log[[1]]

# read daner file
daner_gz <- snakemake@input[[1]]
daner <- read_table2(daner_gz, col_types=cols("SNP"=col_character()))


# remove problem rows
daner_patch <- daner[-problems(daner)$row,]

# check for duplicate SNPs and CPIDs
daner_dups <-
daner_patch %>%
add_count(SNP, name='snp_count') %>%
add_count(CHR, BP, name='cpid_count')

daner_rp <-
daner_dups %>%
filter(snp_count == 1 & cpid_count ==1) %>%
select(-snp_count, -cpid_count)

daner_out_gz <- snakemake@output[[1]]
write_tsv(daner_rp, daner_out_gz)

log_info <- paste("Ricopili check:", "remove", nrow(daner)-nrow(daner_patch), "incomplete rows and", nrow(daner_dups) - nrow(daner_rp), "duplicate markers or positions.")

cat(log_info, file=log_path)