# Convert QC'd Neff daner to GCTA MA format for COJO
# Update value of SE to exactly match log(OR) and P-value

library(readr)
library(dplyr)

daner <- read_tsv(snakemake@input[[1]])


frq_u_col <- names(select(daner, starts_with('FRQ_U')))

# print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"
# calculate se such that 2*pnorm(b/se, lower.tail=F) == P
ma <- daner %>%
mutate(b=log(OR), Z=qnorm(P/2, lower.tail=F)) %>%
transmute(SNP, A1, A2, freq=.data[[frq_u_col]], b, se=b/Z, p=P, N=Neff)

write_delim(ma, snakemake@output[[1]])