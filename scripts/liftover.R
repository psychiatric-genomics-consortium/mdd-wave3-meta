library(dplyr)
library(readr)
library(rtracklayer)


# input daner in hgFROM
daner_hg_from_gz <- snakemake@input$daner
daner_hg_from <- read_table2(daner_hg_from_gz)

# chain file
chain_path <- snakemake@input$chain

# output daner as hgTO
daner_hg_to_gz <- snakemake@output[[1]]

# chain file for hgFROM to hgTO
hgToHg <- import.chain(chain_path)

# create genomic ranges for sumstats
# format chromosomes as 'chrN'
# convert non-autosome numeric chromosomes to text
sumstats_gr_hg_from <- 
with(daner_hg_from,
GRanges(
	seqnames=paste0('chr', case_when(CHR %in% 1:22 ~ as.character(CHR),
                                   CHR == 23 ~ 'X',
                                   CHR == 24 ~ 'Y',
                                   CHR == 25 ~ 'XY',
                                   CHR == 26 ~ 'MT')),
	ranges=IRanges(BP, width=1),
	SNP=SNP
))
	
# liftover from hgFROM to hgTO
sumstats_gr_hg_to <- liftOver(sumstats_gr_hg_from, hgToHg)

# update BP with lifted-over ranges
daner_hg_to <-
daner_hg_from %>%
inner_join(as.data.frame(sumstats_gr_hg_to), by='SNP') %>%
select(CHR, SNP, BP=start, A1, A2, starts_with('FRQ'), INFO, OR, SE, P)

write_tsv(daner_hg_to, path=daner_hg_to_gz)