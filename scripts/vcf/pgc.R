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
medrxiv <- cff$references[[3]]$doi
pmid <- cff$references[[2]]$pmid
pmcid <- cff$references[[2]]$pmcid
doi <- cff$doi
sumstats_url <- cff$references[[1]]$url
code_url <- cff$`repository-code`

analysis_version <- snakemake@wildcards$analysis
cohorts <- snakemake@wildcards$cohorts

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
analysed_cohorts <- read_tsv(snakemake@input$cohorts) |>
  mutate(Neff = 4 * N_cases * N_controls / (N_cases + N_controls))

cat(str_glue("Making sumstats file for cohorts: {snakemake@wildcards$cohort}\n"))

# cohort lists

cohort_list <- str_c(analysed_cohorts$cohort, collapse = ";")
releases_list <- str_c(analysed_cohorts$release, collapse = ";")
studies_list <- str_c(analysed_cohorts$study_name, collapse = "; ")
descriptions_list <- str_c(str_remove(analysed_cohorts$study_description, '"'), collapse = "; ")
reference_list <- str_c(analysed_cohorts$ancestry, collapse = ";")
ncohort <- nrow(analysed_cohorts)

# number of cases and controls
cases_by_cohort <- str_c(analysed_cohorts$N_cases, collapse = ";")
controls_by_cohort <- str_c(analysed_cohorts$N_controls, collapse = ";")
neff_by_cohort <- str_c(round(analysed_cohorts$Neff), collapse = ";")
neff <- round(sum(analysed_cohorts$Neff))
trios_by_cohort <- str_c(rep(0, ncohort), collapse = ";")
processed_by_core <- str_c(analysed_cohorts$data_type == "genotype", collapse=";")

# reference population
ancestries <- snakemake@wildcards$ancestries
reference <- str_to_upper(ancestries)

# number of variants
variants <- nrow(daner)

# header template
cat("Preparing header\n")

header_glues <- readLines(snakemake@input$header_template)
headers <- sapply(header_glues, str_glue)
header <- str_c(headers, collapse = "\n")

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
       INFO, ngt, Neff, Nca, Nco, HetISqt, HetDf, HetPVa) %>%
transmute(`#CHROM`=CHR, POS=BP, ID=SNP, EA=A1, NEA=A2,
          BETA = log(OR), SE, PVAL = P, NGT = ngt, FCAS = FRQ_A, FCON = FRQ_U,
          IMPINFO = INFO, NEFF = Neff,
          NCAS = Nca, NCON = Nco,
          HETI = HetISqt, HETDF = HetDf, HETPVAL = HetPVa)

out <- snakemake@output[[1]]

cat(str_glue("Writing sumstats to {out}"))

cat(header, "\n", file = out)

write_tsv(pgc_sumstats, file = out, append = TRUE, col_names = TRUE)
