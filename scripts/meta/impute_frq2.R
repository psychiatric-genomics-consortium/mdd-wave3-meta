library(dplyr)
library(readr)
library(rtracklayer)
library(stringr)

# log
log_path <- snakemake@log[[1]]

# read imputation reference files
reference_info <- snakemake@input$ref
reference_dir <- readLines(reference_info)[[1]]

# AFREQ file
pop_afreq <- snakemake@input$afreq

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
novel_cups <- read_table(cup_path, col_names=c('chr', 'start', 'end')) %>%
mutate(chr=str_replace(chr, 'chr', '')) %>%
mutate(chr=case_when(chr == 'X' ~ '23',
                     chr == 'Y' ~ '24',
                     TRUE ~ chr))

# list reference files for given ancestries group
frq2_gz <- paste('*', pop, 'frq2.gz', sep='.')
impute_frq2_files <- list.files(reference_dir, pattern=frq2_gz, full.names=T)

# list reference files for chromosome X for given ancestries group
impute_X_frq2_files <- list.files(file.path(reference_dir, 'chr23'), pattern=frq2_gz, full.names=T)

cat('Loading imputation reference\n')

# read in, merge, and remove markers with duplicate positions
impute_frq2 <-
bind_rows(
lapply(c(impute_frq2_files, impute_X_frq2_files),
	   function(frq2_file)
		 read_table(frq2_file,
					 col_types=cols(SNP = col_character(),
									CHR = col_integer(),
									POS = col_integer(),
									A1 = col_character(),
									A2 = col_character(),
									FA1 = col_double(),
									NCHROBS = col_integer()
))))

cat('Loading other reference\n')

# other reference file
afreq <- read_tsv(pop_afreq) %>%
filter(`#CHROM` %in% as.character(c(1:22, 'X'))) %>%
mutate(CHR=if_else(`#CHROM`=='X', true=23, false=as.numeric(`#CHROM`)))

cat('Making genomic ranges\n')

# make genomic range for reference
impute_gr <- with(impute_frq2, GRanges(seqnames=CHR, ranges=IRanges(POS, width=1), SNP=SNP))

# find SNPs with the same range
duplicate_hits <- as_tibble(findOverlaps(impute_gr, impute_gr)) %>% filter(queryHits != subjectHits)

# other reference genomic range
afreq_gr <- with(afreq, GRanges(seqnames=CHR, ranges=IRanges(POS, width=1), SNP=ID))

other_duplicate_hits <- as_tibble(findOverlaps(afreq_gr, afreq_gr)) %>% filter(queryHits != subjectHits)

cat('Finding overlaps\n')

# find overlaps with other reference
ref_overlaps <- findOverlaps(impute_gr, afreq_gr)

cat('Merging reference files\n')

# positions that are in other ref file
other_freq_nodups <- afreq %>%
dplyr::slice(-unique(c(ref_overlaps@to, other_duplicate_hits$queryHits))) %>%
filter(!ID %in% impute_frq2$SNP) %>%
transmute(SNP=ID, CHR, POS, A1=REF, A2=ALT, FA1=REF_FREQ, NCHROBS=OBS_CT, source='1KG')

# merge references
impute_frq2_nodups <-
impute_frq2 %>% dplyr::slice(-duplicate_hits$queryHits) %>%
mutate(source='HRC') %>%
bind_rows(other_freq_nodups) %>%
filter(between(FA1, qc_maf, 1-qc_maf)) %>%
arrange(CHR, POS)

cat('Finishing\n')

saveRDS(impute_frq2_nodups, impute_frq2_rds)

log_info <- str_glue("Refinfo: {reference_info}",
"Otherinfo: {pop_afreq}",
"SNPs: {nrow(impute_frq2)}",
"Other SNPs: {nrow(afreq)}",
"Duplicate SNPs: {nrow(duplicate_hits)}",
"Removed SNPs: {nrow(impute_frq2)-length(unique(duplicate_hits$queryHits))}",
"Extra SNPs: {nrow(other_freq_nodups)}",
"Total SNPs: {nrow(impute_frq2_nodups)}",
"Output: {impute_frq2_rds}",
.sep="\n")

cat(log_info)
cat(log_info, file=log_path)