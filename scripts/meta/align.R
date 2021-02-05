library(dplyr)
library(readr)
library(stringr)

# log
log_path <- snakemake@log[[1]]

# read daner file
daner_gz <- snakemake@input$daner
daner <- read_table2(daner_gz, col_types=cols("SNP"=col_character()))

# read imputation reference files
reference_info <- snakemake@input$ref
reference_dir <- readLines(reference_info)[[1]]

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
select(-NCHROBS)

gc()

# get names of FRQ_A and FRQ_U columns
frq_a_col <- names(select(daner, starts_with('FRQ_A')))
frq_u_col <- names(select(daner, starts_with('FRQ_U')))

# merge on chromosome and position
daner_aligned <- 
daner %>%
inner_join(impute_frq2 %>% filter(between(FA1, 0.01, 0.99)),
          by=c('CHR'='CHR', 'BP'='POS'), suffix=c('', '.imp')) %>%
# keep rows where alleles match
filter((A1 == A1.imp & A2 == A2.imp ) | (A1 == A2.imp & A2 == A1.imp)) %>%
# remove rows with missing statistics
filter(!is.na(OR) & !is.na(SE) & !is.na(P)) %>%
# remove rows where allele frequencies == 0 or 1
filter(.data[[frq_a_col]] > 0, .data[[frq_a_col]] < 1,
       .data[[frq_u_col]] > 0, .data[[frq_u_col]] < 1) %>%
# filter on INFO
filter(INFO >= 0.6) %>%
# select imputed SNP name
mutate(SNP=SNP.imp) %>%
# remove duplicate SNPs
group_by(SNP) %>%
mutate(count=n()) %>%
ungroup() %>%
filter(count == 1) %>%
select(-ends_with('.imp'), -FA1, -count) %>%
arrange(CHR, BP) %>%
select(CHR, SNP, BP, A1, A2, starts_with('FRQ_A'), starts_with('FRQ_U'), INFO, OR, SE, P, everything())

aligned_gz <- snakemake@output[[1]]

write_tsv(daner_aligned, aligned_gz)

# write log
N_snps <- nrow(daner)
mean_or <- exp(mean(log(daner_aligned$OR)))
mean_se <- mean(daner_aligned$SE)
max_or <- max(c(daner_aligned$OR, 1/daner_aligned$OR))
max_se <- max(daner_aligned$SE)
N_snps_aligned <- nrow(daner_aligned)
sumstats_cohort <- snakemake@wildcards$cohort
sumstats_ancestries <- snakemake@wildcards$ancestries
sumstats_release <- snakemake@wildcards$release

log_table <- data.frame(cohort=sumstats_cohort, ancestries=sumstats_ancestries, release=sumstats_release,
           N_snps=N_snps, N_snps_aligned=N_snps_aligned, mean_OR=round(mean_or, 4), max_OR=round(max_or, 4), mean_SE=round(mean_se, 4), max_SE=round(max_se, 4))

write_tsv(log_table, file=log_path)

