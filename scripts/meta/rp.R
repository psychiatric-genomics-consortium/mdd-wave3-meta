# Check Ricopili output

library(dplyr)
library(readr)
library(rtracklayer)
library(stringr)

# log
log_path <- snakemake@log[[1]]

# read daner file
daner_gz <- snakemake@input[[1]]
daner <- read_table2(daner_gz, col_types=cols("SNP"=col_character()))

# remove problem rows
if(nrow(problems(daner)) > 0) {
	daner_patch <- daner %>%
		dplyr::slice(-problems(daner)$row) %>%
		arrange(CHR, BP)
} else {
	daner_patch <- daner %>%
		arrange(CHR, BP)
}

# check for duplicate SNPs and CPIDs

# create genomic ranges for sumstats and reference
daner_gr <- with(daner_patch, GRanges(seqnames=CHR, ranges=IRanges(BP, width=1), SNP=SNP))

# find duplicate positions in the sumstats
sumstats_dups_idx <- as_tibble(findOverlaps(daner_gr, daner_gr)) %>% filter(queryHits != subjectHits) %>% pull(queryHits)

# find duplicate SNPs
snps_dups_idx <- which(duplicated(daner_patch$SNP))

dups_idx <- unique(c(sumstats_dups_idx, snps_dups_idx))

if(length(dups_idx) > 0) {
	daner_rp <- daner_patch %>% dplyr::slice(-dups_idx)
} else {
	daner_rp <- daner_patch
	
}

daner_out_gz <- snakemake@output[[1]]
write_tsv(daner_rp, daner_out_gz)

log_info <- paste("Ricopili check:", "remove", nrow(daner)-nrow(daner_patch), "incomplete rows and", nrow(daner_patch) - nrow(daner_rp), "duplicate markers or positions.")

cat(log_info, file=log_path)