library(readr)
library(dplyr)
library(stringr)
library(tidyr)

levey <- read_tsv(snakemake@input$levey, col_types=cols(CHR.BP=col_character()))

giannakopoulou <- read_tsv(snakemake@input$giannakopoulou, col_types=cols('CHR:POS'=col_character())) %>%
separate(`CHR:POS`, into=c('CHR', 'POS'), convert=TRUE)

gwas_catalog <- read_tsv(snakemake@input$catalog) %>%
filter(!is.na(CHR_ID)) %>%
mutate(CHR=if_else(CHR_ID == 'X', true=23, false=as.numeric(CHR_ID)),
       POS=as.numeric(CHR_POS)) %>%
filter(!is.na(POS))

bim <- read_tsv(snakemake@input$bim, col_names=c('CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'))

# list all variants
previous_snps <-
bind_rows(
select(levey, SNP=rsid),
select(giannakopoulou, SNP),
transmute(gwas_catalog, SNP=paste0('rs', SNP_ID_CURRENT)
)) %>%
distinct()

write_tsv(previous_snps, snakemake@output[[1]], col_names=FALSE)