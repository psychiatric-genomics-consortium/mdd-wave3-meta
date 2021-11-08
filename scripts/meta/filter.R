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
daner <- read_table(daner_gz, col_types=cols("SNP"=col_character()))

N_snps <- nrow(daner)
logging(glue("{N_snps} SNPs read in."))

# get ancestries superpopulation
pop <- toupper(snakemake@wildcards$ancestries)

# dentist outliers
dentist_txt <- snakemake@input$dentist
outliers <- read_table(dentist_txt, col_names='SNP')

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

# Calculating checks

daner_check <- daner %>%
# shorter variable names for frequency columns
mutate(frq_a=.data[[frq_a_col]],
	   frq_u=.data[[frq_u_col]])%>%
# determine MAF
mutate(maf_a=if_else(frq_a <= 0.5, true=frq_a, false=1-frq_a),
	   maf_u=if_else(frq_u <= 0.5, true=frq_u, false=1-frq_u)) %>%
# check for extreme allele frequency differences
mutate(M=(frq_u + FA1.ref)/2,
	   S=(frq_u^2 + FA1.ref^2)/2,
	   Fst=(S-M^2)/(M*(1-M)),
	   D=2*2*(S-M^2)/((2-1)*(1+2*S-2*M)))

median_fst <- median(daner_check$Fst, na.rm=T)
max_fst <- max(daner_check$Fst, na.rm=T)
var_fst <- var(daner_check$Fst, na.rm=T)
logging(glue('Fst: median = {signif(median_fst, 3)}, max = {signif(max_fst, 3)}, var = {signif(var_fst, 3)}'))
	   
logging('Permutation test for Fst')
permute_frq <- daner_check %>%
transmute(frq1=sample(frq_u, n()),
          frq2=sample(FA1.ref, n())) %>%
mutate(M=(frq1 + frq2)/2,
       S=(frq1^2 + frq2^2)/2,
       Fst=(S-M^2)/(M*(1-M)),
       D=2*2*(S-M^2)/((2-1)*(1+2*S-2*M))) 	  
	   
qc_fst <- median(permute_frq$Fst, na.rm=T)


logging('Applying filters')
daner_qc <- daner_check %>%
# Apply QC checks for MAC, MAF, frq different from ref, info
mutate(QC=case_when(maf_a == 0 | maf_u  == 0 ~ 'MAF',
                    maf_a < qc_maf & maf_u < qc_maf ~ 'MAF',
                    maf_a*n_cases < qc_mac | maf_u*n_controls <= qc_mac ~ 'MAC',
					Fst > qc_fst ~ 'FST',
					frq_u - FA1.ref > qc_diff ~ 'DIFF',
					INFO < qc_info ~ 'INFO',
                    SNP %in% outliers$SNP ~ 'DENTIST',
					TRUE ~ 'PASS'))

qc_counts <- table(daner_qc$QC)

snps_maf <- coalesce(qc_counts['MAF'], 0)
snps_mac <- coalesce(qc_counts['MAC'], 0)
snps_dentist <- coalesce(qc_counts['DENTIST'], 0)
snps_fst <- coalesce(qc_counts['FST'], 0)
snps_diff <- coalesce(qc_counts['DIFF'], 0)
snps_info <- coalesce(qc_counts['INFO'], 0)
snps_pass <- coalesce(qc_counts['PASS'], 0)

logging(glue('{snps_maf} SNPs removed for MAF < {qc_maf}.'))
logging(glue('{snps_mac} SNPs removed for MAC < {qc_mac}.'))
logging(glue('{snps_fst} SNPs removed for Fst > {signif(qc_fst, 4)}'))
logging(glue('{snps_diff} SNP removed for DIFF > {qc_diff}.'))
logging(glue('{snps_info} SNP removed for INFO < {qc_info}.'))
logging(glue('{snps_dentist} SNPs removed for DENTIST outlier'))
					
daner_filtered <- daner_qc %>%
filter(QC == 'PASS') %>%
select(-ends_with('.ref'), -frq_a, -frq_u, -maf_a, -maf_u, -QC, -M, -S, -Fst, -D) %>%
arrange(CHR, BP) %>%
select(CHR, SNP, BP, A1, A2, starts_with('FRQ_A'), starts_with('FRQ_U'), INFO, OR, SE, P, everything())



N_snps_kept <- nrow(daner_filtered)
logging(glue("{N_snps_kept} SNPs remain after applying QC filters. {N_snps-N_snps_kept} removed."))		

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

# median OR, SE limited to MAF > 0.01
daner_filtered_maf <- daner_qc %>% filter(QC == 'PASS', between(frq_a, 0.01, 0.99), between(frq_u, 0.01, 0.99)) %>% select(OR, SE)
median_or01 <- median(daner_filtered_maf$OR)
median_se01 <- median(daner_filtered_maf$SE)

qc_table <- data.frame(cohort=sumstats_cohort, ancestries=sumstats_ancestries, release=sumstats_release,
		                N_cases=n_cases, N_controls=n_controls,
						median_fst=signif(median_fst, 3), max_fst=signif(max_fst, 3),
						var_fst=signif(var_fst, 3),
						snps_kept=N_snps_kept,
						median_OR=round(median_or, 4), max_OR=round(max_or, 4), median_OR01=round(median_or01, 4),
						median_SE=round(median_se, 4), max_SE=round(max_se, 4), median_SE01=round(median_se01, 4))

logging(glue("Writing qc table to {qc_txt}"))
write_tsv(qc_table, file=qc_txt)

