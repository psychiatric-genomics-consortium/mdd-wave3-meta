library(dplyr)
library(readr)
library(tidyr)
library(stringr)

args <- commandArgs(TRUE)

text_gz=args[1]
daner_gz=args[2]
logfile=args[3]

sumstats <- read_table2(text_gz, n_max=1000, comment='##')


sumstats_cols <- c("BETA", "SE", "PVAL",
				   "NGT", "FCAS", "FCON",
				   "R2", "NEFF", "NCAS", "NCON", "DIRE")

sumstats_stats <-
sumstats %>%
separate(INFO, into=sumstats_cols, sep=';') %>%
select(-DIRE) %>%
mutate_at(sumstats_cols[1:10], ~ sapply(str_split(.x, pattern='='), function(x) as.numeric(x[2])))

Ncases <- max(sumstats_stats$NCAS)
Ncontrols <- max(sumstats_stats$NCON)

FRQ_A_col <- paste0('FRQ_A_', Ncases)
FRQ_U_col <- paste0('FRQ_U_', Ncontrols)

daner <-
sumstats_stats %>%
mutate(OR=exp(BETA)) %>%
select(CHR=`#CHROM`, SNP=ID, BP=POS, A1=ALT, A2=REF,
  	   !!FRQ_A_col:=FCAS, !!FRQ_U_col:=FCON, INFO=R2,
	   OR, SE, P=PVAL, ngt=NGT)
	   
write_tsv(daner, daner_gz)