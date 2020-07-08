library(dplyr)
library(readr)
library(stringr)

# read daner file
daner_gz <- snakemake@input$daner
daner <- read_table2(daner_gz)

# read imputation reference files
reference_info <- snakemake@input$ref
reference_dir <- readLines(reference_info)[[1]]

# list reference files for given ancestries group
# this needs to be made into a parameter
impute_frq2_files <- list.files(reference_dir, '*.EUR.frq2.gz', full.names=T)

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
select(-FA1, -NCHROBS)

gc()

# merge on chromosome and position
daner_aligned <- 
daner %>%
left_join(impute_frq2, by=c('CHR'='CHR', 'BP'='POS'), suffix=c('', '.imp')) %>%
# keep rows where alleles match or where variant is not in the panel
filter(is.na(SNP.imp) | ((A1 == A1.imp & A2 == A2.imp ) | (A1 == A2.imp & A2 == A1.imp))) %>%
# select imputed SNP name, or original if no match
mutate(SNP=coalesce(SNP.imp, SNP)) %>%
select(-ends_with('.imp'))

aligned_gz <- snakemake@output[[1]]

write_tsv(daner_aligned, aligned_gz)
