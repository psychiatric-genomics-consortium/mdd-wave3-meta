library(dplyr)
library(readr)
library(stringr)
library(yaml)

params <- commandArgs(TRUE)

text_gz <- params[1]
daner_gz <- params[2]
logfile <- params[3]

# reference directory
config <- read_yaml('config.yaml')
reference_dir <- config$refdir

# input sumstats
sumstats <- read_table(text_gz)

impute_frq2_files <- list.files(reference_dir, pattern=paste('*', 'EUR', 'frq2.gz', sep='.'), full.names=T)

impute_frq2 <-
bind_rows(
lapply(impute_frq2_files,
     function(frq2_file)
     read_table(frq2_file,
           col_types=cols(SNP = col_character(),
                  CHR = col_integer(),
                  POS = col_integer(),
                  A1 = col_character(),
                  A2 = col_character(),
                  FA1 = col_double(),
                  NCHROBS = col_integer()
)))) %>%
select(-NCHROBS)

# get marker IDs from reference panel
sumstats_snps <- sumstats %>%
  mutate(EA = str_to_upper(EA), NEA = str_to_upper(NEA)) %>%
  inner_join(impute_frq2,
             by=c("Chromosome" = "CHR", "Position" = "POS"),
             multiple = "any") %>%
  filter((EA == A1 & NEA == A2) | (EA == A2 & NEA == A1))

n_cases <- sumstats %>%
  summarize(Ncases = max(Ncase)) %>%
  pull(Ncases)

n_controls <- sumstats %>%
  summarize(Ncontrols = max(Ncontrol)) %>%
  pull(Ncontrols)

frq_a_col <- str_c("FRQ", "A", n_cases, sep = "_")
frq_u_col <- str_c("FRQ", "U", n_controls, sep = "_")

daner <- sumstats_snps %>%
         mutate(OR = exp(BETA), INFO = 1, Neff_half = Neff / 2) %>%
  select(CHR = Chromosome, SNP, BP = Position,
         A1 = EA, A2 = NEA,
         !!frq_a_col := EAF, !!frq_u_col := EAF, INFO,
         OR, SE, P, Nca = Ncase, Nco = Ncontrol, Neff_half)
         
write_tsv(daner, daner_gz)