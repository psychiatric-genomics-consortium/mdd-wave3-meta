library(dplyr)
library(readr)
library(rtracklayer)
library(stringr)
library(glue)

# log
log_path <- snakemake@log[[1]]
logging <- function(msg, append=TRUE) {
	cat(msg, "\n")
	cat(msg, "\n", file=log_path, append=append)
}

# read daner file
daner_gz <- snakemake@input$daner
logging(glue("Reading in sumstats {daner_gz}"), append=FALSE)
daner <- read_table2(daner_gz, col_types=cols("SNP"=col_character()))

N_snps <- nrow(daner)
logging(glue("{N_snps} SNPs read in."))


# read imputation reference files
impute_frq2_rds <- snakemake@input$ref
logging(glue("Reading in reference {impute_frq2_rds}"))
impute_frq2 <- readRDS(impute_frq2_rds)
logging(glue("{nrow(impute_frq2)} SNPs in reference file"))

# get ancestries superpopulation
pop <- toupper(snakemake@wildcards$ancestries)

# QC paramaters
qc_maf <- snakemake@params$maf
qc_info <- snakemake@params$info
qc_mac <- snakemake@params$mac
qc_secure <- snakemake@params$secure_frq
qc_diff <- snakemake@params$diff_frq

logging("QC parameters")
logging(glue("MAF: {qc_maf} (minor allele frequency)"))
logging(glue("INFO: {qc_info} (imputation score)"))
logging(glue("MAC: {qc_mac} (minor allele count)"))
logging(glue("Secure: {qc_secure} (frequency around 50% to remove strand ambiguous SNPs)"))
logging(glue("Diff: {qc_diff} (difference in allele frequency with reference)"))

# get names of FRQ_A and FRQ_U columns
frq_a_col <- names(select(daner, starts_with('FRQ_A')))
frq_u_col <- names(select(daner, starts_with('FRQ_U')))

n_cases <- as.numeric(str_split(frq_a_col, pattern='_')[[1]][3])
n_controls <- as.numeric(str_split(frq_u_col, pattern='_')[[1]][3])

# create genomic ranges for sumstats and reference
logging("Determining genomic ranges")
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
	encode <- stringi::stri_replace_all_fixed(d, c('+', '-'), c('P', 'M'), vectorise_all=FALSE)
	decode <- stringi::stri_replace_all_fixed(encode, c('P', 'M'), c('-', '+'), vectorise_all=FALSE)
	return(decode)
}

# Merge and filter
logging('Merging with reference panel')
daner_merged <-
daner %>% dplyr::slice(sumstats_impute_align_idx_nodups$queryHits) %>%
bind_cols(impute_frq2 %>% dplyr::slice(sumstats_impute_align_idx_nodups$subjectHits) %>% rename_with(~ paste0(., '.imp'))) %>%
# shorter variable names for frequency columns
mutate(frq_a=.data[[frq_a_col]],
	   frq_u=.data[[frq_u_col]])%>%
# flipped allele values
mutate(A1.flip=flip(A1), A2.flip=flip(A2))

N_snps_merged <- nrow(daner_merged)
logging(glue("{N_snps_merged} SNPs after merge with reference. {N_snps-N_snps_merged} removed."))

rm(daner, impute_frq2)
	   


logging("Checking SNP effects")
daner_nomiss <- 
daner_merged %>%	
# remove rows with missing statistics
filter(!is.na(OR) & !is.na(SE) & !is.na(P)) 

N_snps_nomiss <- nrow(daner_nomiss)
logging(glue("{N_snps_merged-N_snps_nomiss} SNPs removed for missing OR, SE, or P."))

rm(daner_merged)

logging("Checking for unambiguous strand flips")
daner_unambiguous_flips <- daner_nomiss %>%
filter(A1 != A2.flip) %>%
filter((A1.flip == A1.imp & A2.flip == A2.imp) |
       (A1.flip == A2.imp & A2.flip == A1.imp))

N_snps_unambiguous_flips <- nrow(daner_unambiguous_flips)
logging(glue("{N_snps_unambiguous_flips} SNPs with unambiguous strand flips detected"))

if(N_snps_unambiguous_flips > 0) {
	logging("Align while resolving strand orientation")
	resolve_strand <- TRUE
} else {
	logging("Align assuming matching strand orientation")
	resolve_strand <- FALSE
}

# recode to reference panel
logging("Aligning to reference panel")
daner_matching <- daner_nomiss %>%
# remove non-matchable alleles
filter((A1 == A1.imp & A2 == A2.imp) |
	   (A1 == A2.imp & A2 == A1.imp) |
	   (A1.flip == A1.imp & A2.flip == A2.imp) |
	   (A1.flip == A2.imp & A2.flip == A1.imp))
	   
N_snps_matching <- nrow(daner_matching)
logging(glue("{N_snps_merged} SNPs match reference alleles. {N_snps_nomiss-N_snps_matching} removed."))

rm(daner_nomiss)

if(resolve_strand) {
	daner_aligned <- daner_matching %>%
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
							 TRUE ~ NA))
} else {
	daner_aligned <- daner_matching %>%
	mutate(to_turn=case_when(
							# alleles don't need to be turned
							(A1 == A1.imp & A2 == A2.imp) |
							(A1.flip == A1.imp & A2.flip == A2.imp) ~ FALSE,
							 # alleles that need to be turned
							(A1 == A2.imp & A2 == A1.imp) |
							(A1.flip == A2.imp & A2.flip == A1.imp) ~ TRUE,
							TRUE ~ NA))
}

N_snps_turned <- length(which(daner_aligned$to_turn))
N_snps_unturned <- length(which(!daner_aligned$to_turn))
N_snps_unresolved <- length(which(is.na(daner_aligned$to_turn)))
							 
logging(glue("{N_snps_unturned} SNPs with matching effect direction. {N_snps_turned} SNPs turned to match. {N_snps_unresolved} SNPs with ambiguous alleles/frequencies removed."))		
				 
rm(daner_matching)

daner_turned <- daner_aligned %>%
# flip alleles, effects, frequencies to match reference 
filter(!is.na(to_turn)) %>%
mutate(A1=A1.imp,,
	   A2=A2.imp,
	   OR=if_else(to_turn, true=exp(-log(OR)), false=OR),
	   !!frq_a_col := if_else(to_turn, true=1-frq_a, false=frq_a),
	   !!frq_u_col := if_else(to_turn, true=1-frq_u, false=frq_u)) %>%
mutate(across(matches('Direction'), ~if_else(to_turn, true=flip_direction(.), false=.)))

N_snps_aligned <- nrow(daner_turned)

rm(daner_aligned)

daner_filtered <- daner_turned %>%   
# remove rows with small minor allele counts and frequencies
mutate(maf_a=if_else(frq_a <= 0.5, true=frq_a, false=1-frq_a),
	   maf_u=if_else(frq_u <= 0.5, true=frq_u, false=1-frq_u)) %>%
filter(maf_a*n_cases >= qc_mac & maf_u*n_controls >= qc_mac) %>%
filter(maf_a >= qc_maf | maf_u >= qc_maf) %>%
# remove rows where frequency differs greatly from reference
filter(abs(.data[[frq_u_col]] - FA1.imp) <= qc_diff) %>% 
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
logging(glue("{N_snps_kept} SNPs after applying QC filters. {N_snps_aligned-N_snps_kept} removed."))		
rm(daner_turned)

aligned_gz <- snakemake@output$daner

logging(glue("Writing daner to {aligned_gz}"))
write_tsv(daner_filtered, aligned_gz)

# write qc stats
qc_txt <- snakemake@output$snp_counts

median_or <- median(daner_filtered$OR)
mean_se <- mean(daner_filtered$SE)
max_or <- max(c(daner_filtered$OR, 1/daner_filtered$OR))
max_se <- max(daner_filtered$SE)
sumstats_cohort <- snakemake@wildcards$cohort
sumstats_ancestries <- snakemake@wildcards$ancestries
sumstats_release <- snakemake@wildcards$release

qc_table <- data.frame(cohort=sumstats_cohort, ancestries=sumstats_ancestries, release=sumstats_release,
		                N_cases=n_cases, N_controls=n_controls,
						snps=N_snps,
						snps_merged=N_snps_merged,
						snps_unambiguous_flips=N_snps_unambiguous_flips,
						snps_matching=N_snps_matching,
						snps_turned=N_snps_turned,
						snps_unresolved=N_snps_unresolved,
						snps_aligned=N_snps_aligned,
						snps_kept=N_snps_kept,
						median_or=round(median_or, 4), max_OR=round(max_or, 4),
						mean_SE=round(mean_se, 4), max_SE=round(max_se, 4))

logging(glue("Writing qc table to {qc_txt}"))
write_tsv(qc_table, file=qc_txt)

