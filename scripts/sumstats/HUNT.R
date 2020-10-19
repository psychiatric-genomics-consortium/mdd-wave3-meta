library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(parallel)

args <- commandArgs(TRUE)

text_gz=args[1]
daner_gz=args[2]
logfile=args[3]

cat(paste('Counting SNPs', '\n'))
snp_count <- length(count.fields(text_gz, comment.char='#'))

chunk_size <- 500000

chunk_starts <- seq(0, snp_count, by=chunk_size)

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

parse_info <- function(skip, chunk_size, n_total) {
		
	cat(paste('Opening', text_gz, 'at line', skip+1, '\n'))
	sumstats <- read_table2(text_gz, comment='#', skip=skip, n_max=chunk_size,
				col_types='cicccccc', col_names=c('#CHROM', 'POS','ID','REF','ALT','QUAL','FILTER','INFO'))
	
	sumstats_info <- matrix(NA, nrow=nrow(sumstats), ncol=length(sumstats_nums), dimnames=list(NULL, sumstats_nums))
	
	M <- nrow(sumstats)
	for(i in seq.int(M)) {
		if(i %% 100000 == 1) {
			cat(paste0(i+skip, '/', n_total, '\n'))
		}
		sumstats_info[i,] <- as.numeric(assign_key_values(sumstats$INFO[i])[sumstats_nums])
	}
	
	# merge sumstats with the SNP columns
	# convert stats columns to numeric
	sumstats_stats <-
	sumstats %>%
	select(-QUAL, -FILTER, -INFO) %>%
	bind_cols(as_tibble(sumstats_info)) %>%
	filter(between(FCAS, 0.01, .99) | between(FCON, 0.01, 0.99))
	
	return(sumstats_stats)
}

sumstats <- bind_rows(mclapply(chunk_starts, parse_info, chunk_size=chunk_size, n_total=snp_count, mc.cores=16))

Ncases <- max(sumstats$NCAS)
Ncontrols <- max(sumstats$NCON)

FRQ_A_col <- paste0('FRQ_A_', Ncases)
FRQ_U_col <- paste0('FRQ_U_', Ncontrols)

cat(paste('Formatting daner', '\n'))
daner <-
sumstats %>%
mutate(CHR=if_else(`#CHROM` == 'X', true=23, false=as.numeric(`#CHROM`)),
	   OR=exp(BETA)) %>%
select(CHR, SNP=ID, BP=POS, A1=ALT, A2=REF,
		 !!FRQ_A_col:=FCAS, !!FRQ_U_col:=FCON, INFO=R2,
	   OR, SE, P=PVAL, ngt=NGT)

cat(paste('Writing daner', '\n'))   
write_tsv(daner, daner_gz)