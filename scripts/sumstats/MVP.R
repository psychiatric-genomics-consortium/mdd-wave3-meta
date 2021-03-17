library(dplyr)
library(readr)
library(yaml)

args <- commandArgs(TRUE)

text_gz=args[1]
daner_gz=args[2]
logfile=args[3]

sumstats <- read_table2(text_gz,
						col_types=cols(
						rsid = col_character(),
						`CHR:BP` = col_character(),
						CHR = col_double(),
						BP = col_double(),
						A1 = col_character(),
						A2 = col_character(),
						`log(OR)` = col_double(),
						SE = col_double(),
						P = col_double()
						))

# reference directory
config <- read_yaml('config.yaml')
reference_dir <- config$refdir

impute_frq2_files <- list.files(reference_dir, pattern=paste('*', 'EUR', 'frq2.gz', sep='.'), full.names=T)

impute_frq2 <-
bind_rows(
lapply(impute_frq2_files,
	   function(frq2_file)
		 read_table2(frq2_file,
					 col_types=cols(SNP = col_character(),
									CHR = col_integer(),
									POS = col_integer(),
									A1 = col_character(),
									A2 = col_character(),
									FA1 = col_double(),
									NCHROBS = col_integer()
)))) %>%
select(-NCHROBS)

gc()


sumstats_frq <- sumstats %>%
inner_join(impute_frq2, by=c('CHR', 'BP'='POS'), suffix=c('', '.imp')) %>%
# keep rows where alleles match
filter((A1 == A1.imp & A2 == A2.imp ) | (A1 == A2.imp & A2 == A1.imp))

Ncases <- 83810
Ncontrols <- 166405

FRQ_A_col <- paste0('FRQ_A_', Ncases)
FRQ_U_col <- paste0('FRQ_U_', Ncontrols)

cat(paste('Formatting daner', '\n'))
daner <-
sumstats_frq %>%
transmute(CHR, SNP=rsid, BP, A1, A2,
		 !!FRQ_A_col:=if_else(A1 == A1.imp, true=FA1, false=1-FA1),
		 !!FRQ_U_col:=if_else(A1 == A1.imp, true=FA1, false=1-FA1),
		 INFO=1, OR=exp(`log(OR)`), SE, P)

cat(paste('Writing daner', '\n'))   
write_tsv(daner, daner_gz)