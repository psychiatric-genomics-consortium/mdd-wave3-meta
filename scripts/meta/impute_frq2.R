library(dplyr)
library(readr)
library(stringr)

# log
log_path <- snakemake@log[[1]]

# read imputation reference files
reference_info <- snakemake@input$ref
reference_dir <- readLines(reference_info)[[1]]

# output
impute_frq2_rds <- snakemake@output[[1]]

# QC paramaters
qc_maf <- snakemake@params$maf

# get ancestries superpopulation
pop <- toupper(snakemake@wildcards$ancestries)

# list reference files for given ancestries group
impute_frq2_files <- list.files(reference_dir, pattern=paste('*', pop, 'frq2.gz', sep='.'), full.names=T)

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
select(-NCHROBS) %>%
filter(between(FA1, qc_maf, 1-qc_maf))

saveRDS(impute_frq2, impute_frq2_rds)