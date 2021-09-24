# Create VCF-like PGC sumstats file with header

library(dplyr)
library(readr)
library(stringr)
library(yaml)
library(lubridate)

filedate <- format(now(), "%Y-%Om-%d")

# read analyst from config file
analyst <- snakemake@config$analyst

# manuscript and DOI information
biorxiv <- ""
pmid <- ""
pmcid <- ""
doi <- ""
url <- ""

version <- snakemake@wildcards$version

methods <- ""
acknowledgments <- ""
abstract <- ""

# Genome build contig
fasta_fai <- read_table2(snakemake@input$fasta_fai, col_names=c('NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH' ))

# open daner file
daner <- read_tsv(snakemake@input$daner, col_types=cols(SNP=col_character()))

# format contig
daner_chr <- daner %>% select(CHR) %>% distinct(CHR) %>% mutate(CHR=as.character(CHR)) %>% pull(CHR)


contig <- paste(fasta_fai %>%  filter(NAME %in% daner_chr) %>% mutate(contig=str_glue("##contig=<ID={NAME},length={LENGTH}>")) %>% pull(contig), collapse='\n')

# get number of cases and controls from the header
frq_cols_split <- str_split(str_subset(names(daner), 'FRQ'), '_')
ncase = as.numeric(last(frq_cols_split[[1]]))
ncontrol = as.numeric(last(frq_cols_split[[2]]))
ntrio <- 0

# open cohort files

# Genotyped cohorts
genotype_cohort_counts <- read_table2(snakemake@input$genotype_cohorts) %>% filter(Dataset != 'SUM')
genotype_cohorts <- str_match(genotype_cohort_counts$Dataset, 'mdd_(.+)_eur')[,2]
genotype_ncases <- genotype_cohort_counts$N_cases
genotype_ncontrols <- genotype_cohort_counts$N_controls
genotype_nsnps <- genotype_cohort_counts$`N-SNPs`

# Sumstats cohorts
sumstats_cohort_counts <- read_table2(snakemake@input$sumstats_cohorts) %>% filter(Dataset != 'SUM')

# remove PGC MDD from sumstats cohorts list
# extract cohort name
# sum sample sizes over cohorts with multiple samples
sumstats_cohort_noPGC_counts <- sumstats_cohort_counts %>% filter(!str_detect(Dataset, 'mdd_MDD\\d{2}')) %>%
mutate(cohort=str_match(Dataset, 'mdd_(.+)\\.eur')[,2]) %>%
group_by(cohort) %>%
summarize(N_cases=sum(N_cases), N_controls=sum(N_controls), `N-SNPs`=max(`N-SNPs`))

sumstats_cohorts <- sumstats_cohort_noPGC_counts$cohort
sumstats_ncases <- sumstats_cohort_noPGC_counts$N_cases
sumstats_ncontrols <- sumstats_cohort_noPGC_counts$N_controls
sumstats_nsnps <- sumstats_cohort_noPGC_counts$`N-SNPs`

# cohort analysed
cohorts <- snakemake@wildcards$cohort


# cohort lists

cohort_list <- paste(c(genotype_cohorts, sumstats_cohorts), collapse=',')
ncohort <- length(c(genotype_cohorts, sumstats_cohorts))

# number of cases and controls
cases_by_cohort <- paste(c(genotype_ncases, sumstats_ncases), collapse=',')
controls_by_cohort <- paste(c(genotype_ncontrols, sumstats_ncontrols), collapse=',')
trios_by_cohort <- rep(0, ncohort)
snps_by_cohort <- paste(c(genotype_nsnps, sumstats_nsnps), collapse=',')
processed_by_core <- rep(c(TRUE, FALSE), c(length(genotype_cohorts), length(sumstats_cohorts)))

# reference population
ancestries <- toupper(snakemake@wildcards$ancestries)

# number of variants
variants <- nrow(daner)

# header template
header_glue <- paste(readLines(snakemake@input$header_template), collapse='\n')

header <- str_glue(header_glue)

