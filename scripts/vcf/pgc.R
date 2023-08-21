# Create VCF-like PGC sumstats file with header

library(dplyr)
library(readr)
library(stringr)
library(yaml)
library(lubridate)
library(readxl)

cat("Formatting sumstats to pgc.gz\n")

filedate <- format(now(), "%Y-%Om-%d")
year <- snakemake@wildcards$year

# read analyst from config file
analyst <- snakemake@config$analyst

# read citation information
cff <- read_yaml(snakemake@input$cff)

# manuscript and DOI information
biorxiv <- cff$references[[3]]$doi
pmid <- cff$references[[2]]$pmid
pmcid <- cff$references[[2]]$pmcid
doi <- cff$doi
sumstats_url <- cff$references[[1]]$url
code_url <- cff$`repository-code`

analysis_version <- snakemake@params$analysis

methods <- "Summary statistics have been QC'd to Neff >= 80% of max(Neff) [Neff = effective sample size]"
acknowledgments <- "The PGC has received funding from the US National Institute of Mental Health (5 U01MH109528-04). Statistical analyses were carried out on the Genetic Cluster Computer (http://www.geneticcluster.org) hosted by SURFsara and financially supported by the Netherlands Scientific Organization (NWO 480-05-003) along with a supplement from the Dutch Brain Foundation and the VU University Amsterdam."
abstract <- cff$abstract

# Genome build contig
cat(str_glue("Reading contig {snakemake@input$fasta_fai}\n"))
fasta_fai <-
    read_table(snakemake@input$fasta_fai,
        col_names = c("NAME", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"))

# open daner file
cat(str_glue("Reading daner {snakemake@input$daner}\n"))
daner <- read_tsv(snakemake@input$daner,
    col_types = cols(SNP = col_character()))

# format contig
daner_chr <- daner %>%
    select(CHR) %>%
    distinct(CHR) %>%
    mutate(CHR = as.character(CHR)) %>%
    pull(CHR)

# map fai "X" -> "23"
contig_lengths <- fasta_fai %>%
mutate(CHR = case_when(NAME == "X" ~ "23",
                      NAME == "XY" ~ "24",
                      NAME == "MT" ~ "25",
                      TRUE ~ NAME)) %>%
filter(CHR %in% daner_chr) %>%
mutate(contig = str_glue("##contig=<ID={CHR},length={LENGTH}>"))

contig <- paste(pull(contig_lengths, contig), collapse="\n")

# get number of cases and controls from the header
frq_cols_split <- str_split(str_subset(names(daner), 'FRQ'), '_')
ncase = as.numeric(last(frq_cols_split[[1]]))
ncontrol = as.numeric(last(frq_cols_split[[2]]))
ntrio <- 0

# open cohort files

# Analysed cohorts
analysis_counts <- read_excel(snakemake@input$basic) %>%
    mutate(cohort = str_match(Dataset, "mdd_(.+)(\\.|_)eur")[, 2])

# check if genotyped cohorts are included as their own sumstats
includes_pgc_mdd <- "MDD49" %in% analysis_counts$cohort

# Genotyped cohorts
# Subset to those analysed
genotype_cohort_counts <-
    read_table(snakemake@input$genotype_cohorts) %>%
    filter(Dataset != "SUM") %>%
    mutate(cohort=str_match(Dataset, 'mdd_(.+)_eur')[,2]) %>%
    filter(includes_pgc_mdd | cohort %in% analysis_counts$cohort)

genotype_cohorts <- genotype_cohort_counts$cohort
genotype_ncases <- genotype_cohort_counts$N_cases
genotype_ncontrols <- genotype_cohort_counts$N_controls
genotype_neffhalf <- genotype_cohort_counts$N_eff_half
genotype_nsnps <- genotype_cohort_counts$`N-SNPs`

genotype_versions <-
    str_match(genotype_cohort_counts$Dataset,
              "_([:alpha:]+-qc[:digit:]*)")[, 2]

# Sumstats cohorts
# Subset to those analyzed
# remove PGC MDD from sumstats cohorts list
sumstats_cohort_counts <-
    read_table(snakemake@input$sumstats_cohorts) %>%
    filter(Dataset != "SUM") %>%
    mutate(cohort = str_match(Dataset, "mdd_([:alnum:]+)\\.eur")[, 2],
           release = str_match(Dataset, "\\.([[:alnum:]_]+)$")[, 2]) %>%
    filter(cohort %in% analysis_counts$cohort) %>%
    filter(!str_detect(Dataset, "mdd_MDD\\d{2}"))

# extract cohort name
# sum sample sizes over cohorts with multiple samples
sumstats_grouped_counts <- sumstats_cohort_counts %>%
group_by(cohort) %>%
mutate(subcohort = row_number()) %>%
mutate(maxN = max(subcohort)) %>%
ungroup() %>%
mutate(cohortN =
    if_else(maxN == 1,
            true = cohort,
            false = paste0(cohort, subcohort)))

sumstats_cohorts <- sumstats_grouped_counts$cohortN
sumstats_versions <- sumstats_grouped_counts$release
sumstats_ncases <- sumstats_grouped_counts$N_cases
sumstats_ncontrols <- sumstats_grouped_counts$N_controls
sumstats_neffhalf <- sumstats_grouped_counts$N_eff_half
sumstats_nsnps <- sumstats_grouped_counts$`N-SNPs`


# cohort analysed
cohorts <- snakemake@wildcards$cohort

cat(str_glue("Making sumstats file for cohorts: {cohorts}\n"))

# cohort lists

cohort_list <- paste(c(genotype_cohorts, sumstats_cohorts), collapse = ",")
versions_list <- paste(c(genotype_versions, sumstats_versions), collapse = ",")
ncohort <- length(c(genotype_cohorts, sumstats_cohorts))

# number of cases and controls
cases_by_cohort <- paste(c(genotype_ncases, sumstats_ncases), collapse = ",")
controls_by_cohort <- paste(c(genotype_ncontrols, sumstats_ncontrols), collapse = ",")
neff_by_cohort <- paste(2*c(genotype_neffhalf, sumstats_neffhalf), collapse=",")
neff <- 2*sum(c(genotype_neffhalf, sumstats_neffhalf))
trios_by_cohort <- paste(rep(0, ncohort), collapse = ",")
snps_by_cohort <- paste(c(genotype_nsnps, sumstats_nsnps), collapse = ",")
processed_by_core <-
    paste(rep(c(TRUE, FALSE),
              c(length(genotype_cohorts),
              length(sumstats_cohorts))), collapse = ",")

# reference population
ancestries <- toupper(snakemake@wildcards$ancestries)

# number of variants
variants <- nrow(daner)

# header template
cat("Preparing header\n")

header_glues <- readLines(snakemake@input$header_template)
headers <- sapply(header_glues, str_glue)
header <- paste(headers, collapse = "\n")

cat("#########################################################\n")
cat("#########################################################\n")
cat("#########################################################\n")
cat("\n\n\n")
cat(header)

cat("\n\n\n")
cat("#########################################################\n")
cat("#########################################################\n")
cat("#########################################################\n")

# format table
pgc_sumstats <- daner %>%
select(CHR, BP, SNP, A1, A2, OR, SE,
       FRQ_A = starts_with("FRQ_A"), FRQ_U = starts_with("FRQ_U"), P,
       INFO, ngt, Neff, Nca, Nco, Direction, HetISqt, HetDf, HetPVa) %>%
transmute(`#CHROM`=CHR, POS=BP, ID=SNP, A1=A1, A2=A2,
          BETA = log(OR), SE, PVAL = P, NGT = ngt, FCAS = FRQ_A, FCON = FRQ_U,
          IMPINFO = INFO, NEFF = Neff,
          NCAS = Nca, NCON = Nco, DIRE = Direction,
          HETI = HetISqt, HETDF = HetDf, HETPVAL = HetPVa)

out <- snakemake@output[[1]]

cat(str_glue("Writing sumstats to {out}"))

cat(header, "\n", file = out)

write_tsv(pgc_sumstats, file = out, append = TRUE, col_names = TRUE)
