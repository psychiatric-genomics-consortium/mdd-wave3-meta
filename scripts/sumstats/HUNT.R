library(dplyr)
library(readr)
library(tidyr)
library(stringr)

args <- commandArgs(TRUE)

text_gz=args[1]
daner_gz=args[2]
logfile=args[3]

cat(paste('Opening', text_gz, '\n'))
sumstats <- read_table2(text_gz, comment='##', col_types='cicccccc')

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
cat(paste('Parsing INFO', '\n'))
# this could be done functionally but that requries a lot of intermediate
# memory usage
# sumstats_info <- bind_rows(lapply(sumstats$INFO, assign_key_values))
# pre-allocate the final matrix and use a loop instead

sumstats_info <- matrix(NA, nrow=nrow(sumstats), ncol=length(sumstats_nums), dimnames=list(NULL, sumstats_nums))

M <- nrow(sumstats)
for(i in seq.int(M)) {
	sumstats_info[i,] <- as.numeric(assign_key_values(sumstats$INFO[i])[sumstats_nums])
	if(i %% 100000 == 1) {
		cat(paste0(i, '/', M, '\n'))
	}
}

# merge sumstats with the SNP columns
# convert stats columns to numeric
cat(paste('Merging SNP info', '\n'))
sumstats_stats <-
sumstats %>%
select(-QUAL, -FILTER, -INFO) %>%
bind_cols(sumstats_info)

rm(sumstats); rm(sumstats_info); gc()

Ncases <- max(sumstats_stats$NCAS)
Ncontrols <- max(sumstats_stats$NCON)

FRQ_A_col <- paste0('FRQ_A_', Ncases)
FRQ_U_col <- paste0('FRQ_U_', Ncontrols)

cat(paste('Formatting daner', '\n'))
daner <-
sumstats_stats %>%
mutate(CHR=if_else(`#CHROM` == 'X', true=23, false=as.numeric(`#CHROM`)),
       OR=exp(BETA)) %>%
filter(between(FCAS, 0.01, 0.99) & between(FCON, 0.01, 0.99)) %>%
select(CHR, SNP=ID, BP=POS, A1=ALT, A2=REF,
  	   !!FRQ_A_col:=FCAS, !!FRQ_U_col:=FCON, INFO=R2,
	   OR, SE, P=PVAL, ngt=NGT)

cat(paste('Writing daner', '\n'))   
write_tsv(daner, daner_gz)