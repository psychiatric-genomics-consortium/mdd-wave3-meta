library(dplyr)
library(readr)
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
logging(glue("DIFF: {qc_diff} (difference in allele frequency with reference)"))

# get names of FRQ_A and FRQ_U columns
frq_a_col <- names(select(daner, starts_with('FRQ_A')))
frq_u_col <- names(select(daner, starts_with('FRQ_U')))

n_cases <- as.numeric(str_split(frq_a_col, pattern='_')[[1]][3])
n_controls <- as.numeric(str_split(frq_u_col, pattern='_')[[1]][3])

# Filter
logging('Applying filters')
daner_qc <- daner %>%
# shorter variable names for frequency columns
mutate(frq_a=.data[[frq_a_col]],
	   frq_u=.data[[frq_u_col]])%>%
# remove rows with small minor allele counts and frequencies
mutate(maf_a=if_else(frq_a <= 0.5, true=frq_a, false=1-frq_a),
	   maf_u=if_else(frq_u <= 0.5, true=frq_u, false=1-frq_u)) %>%
# Apply QC checks for MAC, MAF, frq different from ref, info
mutate(QC=case_when(maf_a == 0 | maf_u  == 0 ~ 'MAF',
                    maf_a < qc_maf & maf_u < qc_maf ~ 'MAF',
                    maf_a*n_cases < qc_mac | maf_u*n_controls <= qc_mac ~ 'MAC',
					frq_u - FA1.ref > qc_diff ~ 'DIFF',
					INFO < qc_info ~ 'INFO',
					TRUE ~ 'PASS'))

qc_counts <- table(daner_qc$QC)

snps_maf <- coalesce(qc_counts['MAF'], 0)
snps_mac <- coalesce(qc_counts['MAC'], 0)
snps_diff <- coalesce(qc_counts['DIFF'], 0)
snps_info <- coalesce(qc_counts['INFO'], 0)
snps_pass <- coalesce(qc_counts['PASS'], 0)

logging(glue('{snps_maf} SNP removed for MAF.'))
logging(glue('{snps_mac} SNP removed for MAC.'))
logging(glue('{snps_diff} SNP removed for DIFF.'))
logging(glue('{snps_info} SNP removed for INFO.'))
					
daner_filtered <- daner_qc %>%
filter(QC == 'PASS') %>%
select(-ends_with('.ref'), -frq_a, -frq_u, -maf_a, -maf_u, -QC) %>%
arrange(CHR, BP) %>%
select(CHR, SNP, BP, A1, A2, starts_with('FRQ_A'), starts_with('FRQ_U'), INFO, OR, SE, P, everything())



N_snps_kept <- nrow(daner_filtered)
logging(glue("{N_snps_kept} SNPs after applying QC filters. {N_snps-N_snps_kept} removed."))		

filtered_gz <- snakemake@output$daner

logging(glue("Writing daner to {filtered_gz}"))
write_tsv(daner_filtered, filtered_gz)

# write qc stats
qc_txt <- snakemake@output$snp_counts

median_or <- median(daner_filtered$OR)
median_se <- median(daner_filtered$SE)
max_or <- max(c(daner_filtered$OR, 1/daner_filtered$OR))
max_se <- max(daner_filtered$SE)
sumstats_cohort <- snakemake@wildcards$cohort
sumstats_ancestries <- snakemake@wildcards$ancestries
sumstats_release <- snakemake@wildcards$release

qc_table <- data.frame(cohort=sumstats_cohort, ancestries=sumstats_ancestries, release=sumstats_release,
		                N_cases=n_cases, N_controls=n_controls,
						snps_kept=N_snps_kept,
						median_or=round(median_or, 4), max_OR=round(max_or, 4),
						median_SE=round(median_se, 4), max_SE=round(max_se, 4))

logging(glue("Writing qc table to {qc_txt}"))
write_tsv(qc_table, file=qc_txt)

