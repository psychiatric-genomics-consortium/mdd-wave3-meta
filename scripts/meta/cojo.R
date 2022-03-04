# Merge COJO regions together and compare to clumped results,
# regions plots to create list of independent SNPs

library(dplyr)
library(readr)
library(stringr)
library(glue)
library(rtracklayer)

# log
log_path <- snakemake@log[[1]]
logging <- function(msg, append=TRUE) {
    cat(msg, "\n")
    cat(msg, "\n", file=log_path, append=append)
}

regions_jma_cojo <- snakemake@input$cojo

daner_gz <- snakemake@input$daner

clump_txt <- snakemake@input$clump

regions_bim <- snakemake@input$bim

logging("Cojo analysis", append=FALSE)
logging(glue("Sumstats: {daner_gz}"))
logging(glue("Clump file: {clump_txt}"))
logging(glue("COJO regions: {length(regions_jma_cojo)}"))


# Read in COJO jama files
read_cojo <- function(cojo_file) {
    # parse region from filename
    regions_matched <- str_match(basename(cojo_file), "([[:digit:]]+):([[:digit:]]+)-([[:digit:]]+)")
    # read in the file
    cojo_stats <- read_tsv(cojo_file, col_types=cols('SNP'=col_character(), 'refA'=col_character())) %>%
    mutate(snp_idx=row_number(),
           range.left=as.numeric(regions_matched[,3]),
           range.right=as.numeric(regions_matched[,4]))
    return(cojo_stats)
}

cojo_list <- lapply(regions_jma_cojo, read_cojo)
cojo <- bind_rows(cojo_list) 

# bim files from LD reference
bims_list <- lapply(regions_bim, read_tsv, col_names=c('CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'), col_types=cols('SNP'=col_character()))
bims <- bind_rows(bims_list) %>% arrange(CHR, BP)

# Ricopili clumped loci
clump <- read_table2(clump_txt)

# Daner sumstats
daner <- read_tsv(daner_gz, col_types=cols(SNP=col_character()))

# clumped at <= 5e-8
clumped_gw <- clump %>% filter(P <= 5e-8)

logging(glue("Clumped SNPs: {nrow(clumped_gw)}"))
logging(glue("COJO Selected SNPs: {nrow(cojo)}"))

# clumps with top SNP not in the LD reference (not in bim files)
ref_miss_clump <- clumped_gw %>% filter(!SNP %in% bims$SNP)

# look for regions with singleton SNPs
daner_gw <- daner %>% filter(P <= 5e-8)


# Clumped regions that were not analysed in COJO
# Make genomic ranges
cojo_regions_gr <-  with(cojo, GRanges(seqnames=Chr, ranges=IRanges(start=range.left, end=range.right), SNP=SNP))

clump_regions_gr <- with(clump, GRanges(seqnames=CHR, ranges=IRanges(start=range.left, end=range.right), SNP=SNP))

# Find intersection
cojo_clump_overlaps <- findOverlaps(cojo_regions_gr, clump_regions_gr)

# clumps not in intersection
clumps_not_cojoed <- clump %>% dplyr::slice(-unique(cojo_clump_overlaps@to))

# merge with daner information
cojo_daner <-
cojo %>%
left_join(daner, by='SNP') %>%
select(CHR, SNP, BP, A1, A2,
  starts_with('FRQ'), INFO, OR, SE, P, ngt,
  Direction, HetISqt, HetDf, HetPVa, Nca, Nco, Neff_half,
  bJ, bJ_se, pJ, LD_r,
  range.left, range.right)
  
clump_daner <-
clumps_not_cojoed %>%
select(SNP, range.left, range.right) %>%
inner_join(daner, by='SNP') %>%
select(CHR, SNP, BP, A1, A2,
  starts_with('FRQ'), INFO, OR, SE, P, ngt,
  Direction, HetISqt, HetDf, HetPVa, Nca, Nco, Neff_half,
  range.left, range.right)
  
 cojo_clump_daner <-
 bind_rows(cojo_daner, clump_daner)
  
 # count number of GWsig SNPs
 cojo_clump_regions_gr <-  with(cojo_clump_daner, GRanges(seqnames=CHR, ranges=IRanges(start=range.left, end=range.right), SNP=SNP))
 
 daner_gwsig_gr <- with(daner_gw, GRanges(seqnames=CHR, ranges=IRanges(start=BP, width=1), SNP=SNP))

cojo_clump_n_sig_snps <- countOverlaps(cojo_clump_regions_gr, daner_gwsig_gr)

# enumerate regions
cojo_clump_daner_regions <-
cojo_clump_daner %>%
mutate(n_sig_snps=cojo_clump_n_sig_snps)

# non-singleton regions
cojo_nonsingleton_regions <- 
cojo_clump_daner_regions %>%
filter(n_sig_snps > 1) %>%
arrange(CHR, BP) %>%
group_by(CHR, range.left, range.right) %>%
mutate(region=cur_group_id(), snp_idx=row_number()) %>%
ungroup() %>%
select(region, snp_idx, everything())

# singleton regions
cojo_singleton_regions <- 
cojo_clump_daner_regions %>%
filter(n_sig_snps == 1) %>%
arrange(CHR, BP) %>%
group_by(CHR, range.left, range.right) %>%
mutate(region=nrow(cojo_nonsingleton_regions)+cur_group_id(), snp_idx=row_number()) %>%
ungroup() %>%
select(region, snp_idx, everything())


logging(glue("Singleton regions: {nrow(cojo_singleton_regions)}"))

logging(glue("COJO Final SNPs: {nrow(cojo_nonsingleton_regions)}"))
logging(glue("COJO Final SNPs p <= 5e-8, pJ <= 5e-8: {nrow(filter(cojo_nonsingleton_regions, P <= 5e-8, pJ <= 5e-8))}"))
logging(glue("COJO Final SNPs p > 5e-8, pJ <= 5e-8: {nrow(filter(cojo_nonsingleton_regions, P > 5e-8, pJ <= 5e-8))}"))

output_cojo <- snakemake@output$cojo
write_tsv(cojo_nonsingleton_regions, output_cojo)

output_singletons <- snakemake@output$singletons
write_tsv(cojo_singleton_regions, output_singletons)







