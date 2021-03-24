library(dplyr)
library(readr)
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
qc_mac <- snakemake@params$maf

# Read in QC'd list of SNPs from the imputation panel
impute_frq2 <- readRDS(impute_frq2_rds)

# get names of FRQ_A and FRQ_U columns
frq_a_col <- names(select(daner, starts_with('FRQ_A')))
frq_u_col <- names(select(daner, starts_with('FRQ_U')))

n_cases <- as.numeric(str_split(frq_a_col, pattern='_')[[1]][3])
n_controls <- as.numeric(str_split(frq_u_col, pattern='_')[[1]][3])

# merge on chromosome and position
daner_aligned <- 
daner %>%
inner_join(impute_frq2, by=c('CHR'='CHR', 'BP'='POS'), suffix=c('', '.imp')) %>%
# keep rows where alleles match
filter((A1 == A1.imp & A2 == A2.imp ) | (A1 == A2.imp & A2 == A1.imp)) %>%
# remove rows with missing statistics
filter(!is.na(OR) & !is.na(SE) & !is.na(P)) %>%
# remove rows with small minor allele counts
mutate(frq_a=.data[[frq_a_col]],
       frq_u=.data[[frq_u_col]]) %>%
mutate(maf_a=if_else(frq_a <= 0.5, true=frq_a, false=1-frq_a),
       maf_u=if_else(frq_u <= 0.5, true=frq_u, false=1-frq_u)) %>%
filter(maf_a*n_cases >= qc_mac & maf_u*n_controls >= qc_mac) %>%
# filter on INFO
filter(INFO >= qc_info) %>%
# select imputed SNP name
mutate(SNP=SNP.imp) %>%
# remove duplicate SNPs
add_count(SNP, name="snp_count") %>%
filter(snp_count == 1) %>%
# remove duplicate CHR/BP (cpid)
add_count(CHR, BP, name="cpid_count") %>%
filter(cpid_count == 1) %>%
select(-ends_with('.imp'), -FA1, -snp_count, -cpid_count, -frq_a, -frq_u, -maf_a, -maf_u) %>%
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

