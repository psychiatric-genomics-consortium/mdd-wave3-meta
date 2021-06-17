library(dplyr)
library(readr)
library(rtracklayer)
library(stringr)

# log
log_path <- snakemake@log[[1]]

# read daner file
daner_gz <- snakemake@input$daner
daner <- read_table2(daner_gz, col_types=cols("SNP"=col_character()))

# read imputation reference files
impute_frq2_rds <- snakemake@input$ref

# get ancestries superpopulation
pop <- toupper(snakemake@wildcards$ancestries)

# QC paramaters
qc_maf <- snakemake@params$maf
qc_info <- snakemake@params$info
qc_mac <- snakemake@params$mac

# Read in QC'd list of SNPs from the imputation panel
impute_frq2 <- readRDS(impute_frq2_rds)

# get names of FRQ_A and FRQ_U columns
frq_a_col <- names(select(daner, starts_with('FRQ_A')))
frq_u_col <- names(select(daner, starts_with('FRQ_U')))

n_cases <- as.numeric(str_split(frq_a_col, pattern='_')[[1]][3])
n_controls <- as.numeric(str_split(frq_u_col, pattern='_')[[1]][3])

# create genomic ranges for sumstats and reference
sumstats_gr <- with(daner, GRanges(seqnames=CHR, ranges=IRanges(BP, width=1), SNP=SNP))
impute_gr <- with(impute_frq2, GRanges(seqnames=CHR, ranges=IRanges(POS, width=1), SNP=SNP))

# find duplicate positions in the sumstats
sumstats_dups_idx <- as_tibble(findOverlaps(sumstats_gr, sumstats_gr)) %>% filter(queryHits != subjectHits)

# match up sumstats and imputation panel
sumstats_impute_align_idx <- as_tibble(findOverlaps(sumstats_gr, impute_gr))
sumstats_impute_align_idx_nodups <- sumstats_impute_align_idx %>% filter(!queryHits %in% sumstats_dups_idx$queryHits)

# reference panel alleles
impute_alleles <- unique(c(impute_frq2$A1, impute_frq2$A2))

flip <- function(nucleotide) {
	case_when(nucleotide == 'A' ~ 'T',
	          nucleotide == 'T' ~ 'A',
			  nucleotide == 'C' ~ 'G',
			  nucleotide == 'G' ~ 'C',
			  TRUE ~ nucleotide)
}

flip_strand <- function(allele) {
	str_replace_all(allele, pattern='[ATCG]', replacement=flip)
}

# merge using alignment indices
daner_aligned <-
daner %>% dplyr::slice(sumstats_impute_align_idx_nodups$queryHits) %>%
bind_cols(impute_frq2 %>% dplyr::slice(sumstats_impute_align_idx_nodups$subjectHits) %>% rename_with(~ paste0(., '.imp'))) %>%
# remove rows with missing statistics
filter(!is.na(OR) & !is.na(SE) & !is.na(P)) %>%
# remove alleles not found in the reference
filter(A1 %in% impute_alleles & A2 %in% impute_alleles) %>%
# remove non-matchable alleles
filter((A1 == A1.imp & A2 == A2.imp) |
       (A1 == A2.imp & A2 == A1.imp) |
       (A1 == flip(A1.imp) & A2 == flip(A2.imp)) |
	   (A1 == flip(A2.imp) & A2 == flip(A1.imp))) %>%
# remove rows with small minor allele counts and frequencies
mutate(frq_a=.data[[frq_a_col]],
	   frq_u=.data[[frq_u_col]]) %>%
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
select(-ends_with('.imp'), -frq_a, -frq_u, -maf_a, -maf_u) %>%
arrange(CHR, BP) %>%
select(CHR, SNP, BP, A1, A2, starts_with('FRQ_A'), starts_with('FRQ_U'), INFO, OR, SE, P, everything())

aligned_gz <- snakemake@output[[1]]

write_tsv(daner_aligned, aligned_gz)

# write log

N_snps <- nrow(daner)
N_snps_aligned <- nrow(daner_aligned)
median_or <- median(daner_aligned$OR)
mean_se <- mean(daner_aligned$SE)
max_or <- max(c(daner_aligned$OR, 1/daner_aligned$OR))
max_se <- max(daner_aligned$SE)
sumstats_cohort <- snakemake@wildcards$cohort
sumstats_ancestries <- snakemake@wildcards$ancestries
sumstats_release <- snakemake@wildcards$release

log_table <- data.frame(cohort=sumstats_cohort, ancestries=sumstats_ancestries, release=sumstats_release,
		                N_cases=n_cases, N_controls=n_controls, N_snps=N_snps, N_snps_aligned=N_snps_aligned,
						median_or=round(median_or, 4), max_OR=round(max_or, 4),
						mean_SE=round(mean_se, 4), max_SE=round(max_se, 4))

write_tsv(log_table, file=log_path)

