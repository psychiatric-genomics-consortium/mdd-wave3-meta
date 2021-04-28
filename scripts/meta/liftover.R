library(dplyr)
library(readr)
library(rtracklayer)
library(stringr)

# logfile
log_path <- snakemake@log[[1]]

# input daner in hgFROM
daner_hg_from_gz <- snakemake@input$daner
daner_hg_from <- read_table2(daner_hg_from_gz)

# chain file
chain_path <- snakemake@input$chain

# CUP file
cup_path <- snakemake@input$cups

# output daner as hgTO
daner_hg_to_gz <- snakemake@output[[1]]

# chain file for hgFROM to hgTO
hgToHg <- import.chain(chain_path)

# conversion-unstable positions file of ranges to remove
novel_cups <- read_table2(cup_path, col_names=c('chr', 'start', 'end'))

# create genomic ranges for sumstats (assume scaffold is called 'chrN')
# convert non-autosome numeric chromosomes to text
sumstats_gr_hg_from <- 
with(daner_hg_from,
GRanges(
	seqnames=case_when(CHR == 'chr23' ~ 'chrX',
                       CHR == 'chr24' ~ 'chrY',
                       CHR == 'chr25' ~ 'chrXY',
                       CHR == 'chr26' ~ 'chrMT',
					   TRUE ~ CHR),
	ranges=IRanges(BP, width=1),
	SNP=SNP
))

# create genomic ranges for CUPs
cups_gr <- GRanges(seqname=novel_cups$chr, ranges=IRanges(start=novel_cups$start, end=novel_cups$end))
	
# remove CUPs	
sumstats_gr_hg_from_nocups <- sumstats_gr_hg_from[sumstats_gr_hg_from %outside% cups_gr]

# liftover from hgFROM to hgTO
sumstats_gr_hg_to <- liftOver(sumstats_gr_hg_from_nocups, hgToHg)

# update BP with lifted-over ranges
daner_hg_to <-
daner_hg_from %>%
inner_join(as.data.frame(sumstats_gr_hg_to), by='SNP') %>%
mutate(CHR=str_replace(CHR, 'chr', '')) %>%
select(CHR, SNP, BP=start, A1, A2, starts_with('FRQ'), INFO, OR, SE, P)

write_tsv(daner_hg_to, daner_hg_to_gz)

log_info <- str_glue("Sumstats: {daner_hg_from_gz}",
"SNPs: {nrow(daner_hg_from)}",
"Chain: {chain_path}",
"CUPs: {cup_path}",
"CUP regions: {nrow(novel_cups)}",
"Removed SNPs: {nrow(daner_hg_from)-nrow(daner_hg_to)}",
"Kept SNPs: {nrow(daner_hg_to)}",
"Output: {daner_hg_to_gz}",
.sep="\n")

cat(log_info)
cat(log_info, file=log_path)
