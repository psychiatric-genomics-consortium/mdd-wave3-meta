library(dplyr)
library(readr)
library(tidyr)
library(stringr)

args <- commandArgs(TRUE)

text_gz=args[1]
daner_gz=args[2]
logfile=args[3]

sumstats <- read_table2(text_gz, comment='##')

# data that needs to be converted to numeric
sumstats_nums <- c("BETA", "SE", "PVAL",
				   "NGT", "FCAS", "FCON",
				   "R2", "NEFF", "NCAS", "NCON")
				   

# parse out the INFO column in the format key1=value;key2=value;key3=value;...
# return a vector of values named by key
assign_key_values <- function(INFO) {
	keys_values <- str_split(str_split(INFO, pattern=';')[[1]], pattern='=')
	keys <- sapply(keys_values, first)
	values <- sapply(keys_values, last)
	names(values) <- keys
	return(values)
}

# parse key/values and bind into a data.frame
sumstats_info <- bind_rows(lapply(sumstats$INFO, assign_key_values))

# merge sumstats with the SNP columns
# convert stats columns to numeric
sumstats_stats <-
sumstats %>%
select(-QUAL, -FILTER, -INFO) %>%
bind_cols(sumstats_info) %>%
mutate_at(sumstats_nums, as.numeric)

rm(sumstats); rm(sumstats_info); gc()

Ncases <- max(sumstats_stats$NCAS)
Ncontrols <- max(sumstats_stats$NCON)

FRQ_A_col <- paste0('FRQ_A_', Ncases)
FRQ_U_col <- paste0('FRQ_U_', Ncontrols)

daner <-
sumstats_stats %>%
mutate(OR=exp(BETA)) %>%
filter(between(FCAS, 0.01, 0.99) & between(FCON, 0.01, 0.99)) %>%
select(CHR=`#CHROM`, SNP=ID, BP=POS, A1=ALT, A2=REF,
  	   !!FRQ_A_col:=FCAS, !!FRQ_U_col:=FCON, INFO=R2,
	   OR, SE, P=PVAL, ngt=NGT)
	   
write_tsv(daner, daner_gz)