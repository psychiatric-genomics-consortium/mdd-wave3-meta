library(dplyr)
library(readr)
library(rtracklayer)
library(stringr)

# log
log_path <- snakemake@log[[1]]

# read imputation reference files
reference_info <- snakemake@input$ref
reference_dir <- readLines(reference_info)[[1]]

# CUP file
cup_path <- snakemake@input$cups

# output
impute_frq2_rds <- snakemake@output[[1]]

# QC paramaters
qc_maf <- snakemake@params$maf

# get ancestries superpopulation
pop <- toupper(snakemake@wildcards$ancestries)

# conversion-unstable positions file of ranges to remove
# Reference files uses '23' for X chromosome
novel_cups <- read_table2(cup_path, col_names=c('chr', 'start', 'end')) %>%
mutate(chr=str_replace(chr, 'chr', '')) %>%
mutate(chr=case_when(chr == 'X' ~ '23',
                     chr == 'Y' ~ '24',
                     TRUE ~ chr))

# list reference files for given ancestries group
frq2_gz <- paste('*', pop, 'frq2.gz', sep='.')
impute_frq2_files <- list.files(reference_dir, pattern=frq2_gz, full.names=T)

# list reference files for chromosome X for given ancestries group
impute_X_frq2_files <- list.files(file.path(reference_dir, 'chr23'), pattern=frq2_gz, full.names=T)

# read in, merge, and remove markers with duplicate positions
impute_frq2 <-
bind_rows(
lapply(c(impute_frq2_files, impute_X_frq2_files),
	   function(frq2_file)
		 read_table2(frq2_file,
					 col_types=cols(SNP = col_character(),
									CHR = col_integer(),
									POS = col_integer(),
									A1 = col_character(),
									A2 = col_character(),
									FA1 = col_double(),
									NCHROBS = col_integer()
))))

# make genomic range
impute_gr <- with(impute_frq2, GRanges(seqnames=CHR, ranges=IRanges(POS, width=1), SNP=SNP))

# find SNPs with the same range
duplicate_hits <- as_tibble(findOverlaps(impute_gr, impute_gr)) %>% filter(queryHits != subjectHits)

# find overlaps with CUPs
# create genomic ranges for CUPs
cups_gr <- with(novel_cups, GRanges(seqname=chr, ranges=IRanges(start=start, end=end)))

cups_overlaps <- findOverlaps(impute_gr, cups_gr)

impute_frq2_nodups <-
impute_frq2 %>% dplyr::slice(-duplicate_hits$queryHits)

saveRDS(impute_frq2_nodups, impute_frq2_rds)

log_info <- str_glue("Refinfo: {reference_info}",
"SNPs: {nrow(impute_frq2)}",
"CUPs: {cup_path}",
"CUP regions: {nrow(novel_cups)}",
"Duplicate SNPs: {nrow(duplicate_hits)}",
"CUPs SNPs: {length(cups_overlaps@from)}",
"Removed SNPs: {nrow(impute_frq2)-nrow(impute_frq2_nodups)}",
"Kept SNPs: {nrow(impute_frq2_nodups)}",
"Output: {impute_frq2_rds}",
.sep="\n")

cat(log_info)
cat(log_info, file=log_path)