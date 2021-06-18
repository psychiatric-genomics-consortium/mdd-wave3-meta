library(dplyr)
library(readr)
library(rtracklayer)
library(stringr)

# log
log_path <- snakemake@log[[1]]

# read daner file
daner_gz <- snakemake@input$daner
cat("Reading in", daner_gz, "\n")
daner <- read_table2(daner_gz, col_types=cols("SNP"=col_character()))

# read imputation reference files
impute_frq2_rds <- snakemake@input$ref
cat("Reading in", impute_frq2_rds, "\n")
impute_frq2 <- readRDS(impute_frq2_rds)

# get ancestries superpopulation
pop <- toupper(snakemake@wildcards$ancestries)

# QC paramaters
qc_maf <- snakemake@params$maf
qc_info <- snakemake@params$info
qc_mac <- snakemake@params$mac
qc_secure <- snakemake@params$secure_frq
qc_diff <- snakemake@params$diff_frq

# get names of FRQ_A and FRQ_U columns
frq_a_col <- names(select(daner, starts_with('FRQ_A')))
frq_u_col <- names(select(daner, starts_with('FRQ_U')))

n_cases <- as.numeric(str_split(frq_a_col, pattern='_')[[1]][3])
n_controls <- as.numeric(str_split(frq_u_col, pattern='_')[[1]][3])

# create genomic ranges for sumstats and reference
cat("Determining genomic ranges", "\n")
sumstats_gr <- with(daner, GRanges(seqnames=CHR, ranges=IRanges(BP, width=1), SNP=SNP))
impute_gr <- with(impute_frq2, GRanges(seqnames=CHR, ranges=IRanges(POS, width=1), SNP=SNP))

# find duplicate positions in the sumstats
sumstats_dups_idx <- as_tibble(findOverlaps(sumstats_gr, sumstats_gr)) %>% filter(queryHits != subjectHits)

# match up sumstats and imputation panel
sumstats_impute_align_idx <- as_tibble(findOverlaps(sumstats_gr, impute_gr))
sumstats_impute_align_idx_nodups <- sumstats_impute_align_idx %>% filter(!queryHits %in% sumstats_dups_idx$queryHits)

# Flip alleles. handle arbitrary allele length
flip <- function(allele) {
	allele_lower <- tolower(allele)
	stringi::stri_replace_all_fixed(allele_lower, c('a', 'c', 'g', 't'), c('T', 'G', 'C', 'A'), vectorise_all=FALSE)
}

# flip effect direction strings
flip_direction <- function(d) {
	encode <- stringi::stri_replace_all_fixed(d, c('+', '-'), c('P', 'N'), vectorise_all=FALSE)
	decode <- stringi::stri_replace_all_fixed(encode, c('P', 'N'), c('-', '+'), vectorise_all=FALSE)
	return(decode)
}


# Merge and filter
cat('Merging with reference panel', '\n')
daner_merged <-
daner %>% dplyr::slice(sumstats_impute_align_idx_nodups$queryHits) %>%
bind_cols(impute_frq2 %>% dplyr::slice(sumstats_impute_align_idx_nodups$subjectHits) %>% rename_with(~ paste0(., '.imp'))) %>%
# shorter variable names for frequency columns
mutate(frq_a=.data[[frq_a_col]],
	   frq_u=.data[[frq_u_col]]) %>%	
# remove rows with missing statistics
filter(!is.na(OR) & !is.na(SE) & !is.na(P)) %>%
# flipped allele values
mutate(A1.flip=flip(A1), A2.flip=flip(A2))

N_snps <- nrow(daner)
N_snps_merged <- nrow(daner_merged)
rm(daner, impute_frq2)

# recode to reference panel
cat("Aligning to reference panel", "\n")
daner_aligned <- daner_merged %>%
# remove non-matchable alleles
filter((A1 == A1.imp & A2 == A2.imp) |
	   (A1 == A2.imp & A2 == A1.imp) |
	   (A1.flip == A1.imp & A2.flip == A2.imp) |
	   (A1.flip == A2.imp & A2.flip == A1.imp)) %>%
mutate(to_turn=case_when(
						# unambiguous alleles that don't need to be turned
						A1 != A2.flip &
                            ((A1 == A1.imp & A2 == A2.imp) |
							 (A1.flip == A1.imp & A2.flip == A2.imp))
					     ~ FALSE,
						 # unambiguous alleles that need to be turned
						 A1 != A2.flip & 
						 	((A1 == A2.imp & A2 == A1.imp) |
						 	 (A1.flip == A2.imp & A2.flip == A1.imp))
						~ TRUE,
						# ambiguous alleles that don't need to be turned
						 A1 == A2.flip & 
						 	((A1 == A1.imp & A2 == A2.imp) |
						     (A1.flip == A1.imp & A2.flip == A2.imp)) &
						    abs(frq_u - FA1.imp) <= abs((1-frq_u) - FA1.imp) &
							abs(frq_u - 0.5) > qc_secure / 2 
						~ FALSE,
						# ambiguous alleles that do need to be turned
						 A1 == A2.flip & 
						    ((A1 == A2.imp & A2 == A1.imp) | 
							 (A1.flip == A2.imp & A2.flip == A1.imp)) &
							abs(frq_u - FA1.imp) >= abs((1-frq_u) - FA1.imp) &
							abs(frq_u - 0.5) > qc_secure / 2 
						~ TRUE,
						# ambiguous alleles with unresolvable flipping
						 TRUE ~ NA)) %>%
# flip alleles, effects, frequencies to match reference 
filter(!is.na(to_turn)) %>%
mutate(A1=A1.imp,,
	   A2=A2.imp,
	   OR=if_else(to_turn, true=exp(-log(OR)), false=OR),
	   !!frq_a_col := if_else(to_turn, true=1-frq_a, false=frq_a),
	   !!frq_u_col := if_else(to_turn, true=1-frq_u, false=frq_u)) %>%
mutate(across(matches('Direction'), ~if_else(to_turn, true=flip_direction(.), false=.)))

N_snps_aligned <- nrow(daner_aligned)
rm(daner_merged)

cat('Filtering', '\n')
daner_filtered <- daner_aligned %>%   
# remove rows with small minor allele counts and frequencies
mutate(maf_a=if_else(frq_a <= 0.5, true=frq_a, false=1-frq_a),
	   maf_u=if_else(frq_u <= 0.5, true=frq_u, false=1-frq_u)) %>%
filter(maf_a*n_cases >= qc_mac & maf_u*n_controls >= qc_mac) %>%
filter(maf_a >= qc_maf | maf_u >= qc_maf) %>%
# filter on INFO
filter(INFO >= qc_info) %>%
# select imputation reference SNP name
mutate(SNP=SNP.imp) %>%
# remove duplicate SNPs
filter(!duplicated(SNP)) %>%
select(-ends_with('.imp'), -frq_a, -frq_u, -maf_a, -maf_u, -to_turn) %>%
arrange(CHR, BP) %>%
select(CHR, SNP, BP, A1, A2, starts_with('FRQ_A'), starts_with('FRQ_U'), INFO, OR, SE, P, everything())

N_snps_kept <- nrow(daner_filtered)
rm(daner_aligned)

aligned_gz <- snakemake@output[[1]]

cat("Writing daner", aligned_gz, "\n")
write_tsv(daner_filtered, aligned_gz)

# write log

median_or <- median(daner_filtered$OR)
mean_se <- mean(daner_filtered$SE)
max_or <- max(c(daner_filtered$OR, 1/daner_filtered$OR))
max_se <- max(daner_filtered$SE)
sumstats_cohort <- snakemake@wildcards$cohort
sumstats_ancestries <- snakemake@wildcards$ancestries
sumstats_release <- snakemake@wildcards$release

log_table <- data.frame(cohort=sumstats_cohort, ancestries=sumstats_ancestries, release=sumstats_release,
		                N_cases=n_cases, N_controls=n_controls,
						snps=N_snps,
						snps_merged=N_snps_merged,
						snps_aligned=N_snps_aligned,
						snps_kept=N_snps_kept,
						median_or=round(median_or, 4), max_OR=round(max_or, 4),
						mean_SE=round(mean_se, 4), max_SE=round(max_se, 4))

write_tsv(log_table, file=log_path)

