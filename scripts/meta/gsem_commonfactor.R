
library(GenomicSEM)
library(dplyr)
library(readr)

# load the covariance structure

covstruct <- dget(snakemake@input$covstruct)

# load the sumstats for the region requested
sumstats <- read_tsv(snakemake@input$sumstats) %>%
filter(CHR==snakemake@wildcards$CHR,
       between(BP, snakemake@wildcards$START, snakemake@wildcards$STOP)) %>%
as.data.frame()

gwas <- commonfactorGWAS(covstruc=covstruct,
                         SNPs=sumstats,
                         parallel=FALSE)
                         
write_tsv(gwas, snakemake@output[[1]])
