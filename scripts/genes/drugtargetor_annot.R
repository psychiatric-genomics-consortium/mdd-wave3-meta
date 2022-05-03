## DrugTargetor https://drugtargetor.com/ annotations for LDSC

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(rtracklayer)

# Drug Targetor drug/gene association list
targetor <- read_tsv(snakemake@input$drugtargetor)

# GeneMatrix for gene positional information
genematrix <- read_tsv(snakemake@input$geneMatrix)

# bim file for SNPs to use in making LD Scores
bim <- read_table(snakemake@params$bim, col_names=c('CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'))

# annotation window
annot_window <- snakemake@params$window 

# group drugs together using ATC  <https://www.whocc.no/atc/structure_and_principles/> chemical subgroups (4th level)

# ATC codes have the form "A12BC34" and we want the "A12BC" part
# ATC code is tagged within the `atc` column as `ATC:`, along
# with tags for `NAME` and `CID`, separated by "|". Split out
# the first elment split by "|", remove the "ATC:" prefix,
# then extract using a regex (so we are sure we are matching)
# the expected pattern. Join by gene name with the geneMatrix
# to get start and end positions of the gene

targetor_atc4 <-
targetor %>%
select(atc, gene) %>%
mutate(ATC=str_replace(sapply(str_split(atc, '[|]'), dplyr::first),
                       pattern="ATC:",
                       replace="")) %>%
mutate(ATC4=str_extract(ATC, '[A-Z][0-9]{2}[A-Z]{2}')) %>%
distinct(gene, ATC4) %>%
inner_join(genematrix, by=c('gene'='gene_name')) %>%
mutate(chr=str_replace(hg19g0, "chr", "")) %>%
transmute(ensgid, gene, ATC4,
        CHR=case_when(chr == "X" ~ 23,
                      chr == "Y" ~ 24,
                      chr == "M" ~ 25,
                      TRUE ~ as.numeric(chr)),
        START=g1, END=g2) %>%
filter(!is.na(CHR))
        
# Merge into annotation with the bim file. Make genomic ranges to get
# all SNPs within `window` of the start/end of each gene
targetor_gr <- with(targetor_atc4,
    GRanges(seqnames=CHR,
            ranges=IRanges(start=START-annot_window, 
                           end=END+annot_window)))

bim_gr <- with(bim, GRanges(seqnames=CHR, ranges=IRanges(start=BP, width=1)))

# Find overlaps
targetor_bim_overlaps <- findOverlaps(targetor_gr, bim_gr)

# Merge matching rows from both sets
snps_atc <- bind_cols(
dplyr::slice(targetor_atc4 %>% select(-CHR), targetor_bim_overlaps@from),
dplyr::slice(bim, targetor_bim_overlaps@to)
)

# Pivot wider with rows indexed by SNPs and a column for each annotation
snps_annot <- 
snps_atc %>%
transmute(SNP, ATC4, value=1) %>%
distinct() %>%
pivot_wider(id_cols='SNP', names_from='ATC4', values_from='value')

# merge with bim file and fill in NAs with 0
bim_annot <- bim %>%
select(CHR, SNP, BP, CM) %>%
left_join(snps_annot, by='SNP') %>%
mutate(across(-CHR:-CM, .fns=~if_else(is.na(.x), true=0, false=.x)))

write_tsv(bim_annot, snakemake@output[[1]])