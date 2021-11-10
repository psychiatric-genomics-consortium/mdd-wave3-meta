# Convert daner sumstats to GCTA MA format

library(readr)
library(dplyr)
library(stringr)

daner_gz <- snakemake@input[[1]]
cat('Reading', daner_gz, '\n')
daner <- read_tsv(daner_gz)

# get names of FRQ_A and FRQ_U columns
frq_a_col <- names(select(daner, starts_with('FRQ_A')))
frq_u_col <- names(select(daner, starts_with('FRQ_U')))

n_cases <- as.numeric(str_split(frq_a_col, pattern='_')[[1]][3])
n_controls <- as.numeric(str_split(frq_u_col, pattern='_')[[1]][3])

# check if daner file has Nca/Nco rows
if(all(c('Nca', 'Nco') %in% names(daner))) {
    cat('Found Nca, Nco columns', '\n')
    daner_n <- daner
} else {
    cat('Adding Nca, Nco from header', '\n')
    daner_n <- daner %>% mutate(Nca=n_cases, Nco=n_controls)
}

ma <- daner_n %>%
transmute(SNP, A1, A2, freq=.data[[frq_u_col]],
          beta=signif(log(OR), 6), se=SE, p=P,
          N=round(4*INFO*Nca*Nco/(Nca+Nco))) %>%
filter(freq != 0)
          
out_ma <- snakemake@output[[1]]
cat('Writing', out_ma, '\n')
write_delim(ma, out_ma)