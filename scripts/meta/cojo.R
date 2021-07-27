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
    # read in the file
    cojo_stats <- read_tsv(cojo_file, col_types=cols('SNP'=col_character(), 'refA'=col_character())) %>%
    mutate(snp_idx=row_number())
    return(cojo_stats)
}

cojo_list <- lapply(regions_jma_cojo, read_cojo)
cojo <- bind_rows(cojo_list) 

# region coordinates
regions_matched <- str_match(basename(regions_jma_cojo), "([[:digit:]]+):([[:digit:]]+)-([[:digit:]]+)")
regions <- tibble(Chr=as.numeric(regions_matched[,2]),
                  range.left=as.numeric(regions_matched[,3]),
                  range.right=as.numeric(regions_matched[,4])
)

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

# regions with empty COJO files
ref_miss_cojo <- regions %>% dplyr::slice(which(sapply(cojo_list, nrow) ==0))

ref_miss_cojo_clump <- clumped_gw %>% inner_join(ref_miss_cojo, by=c('CHR'='Chr', 'range.left', 'range.right'))

# count number of significant SNPs in each region
region_snp_counts <- plyr::adply(regions, 1, function(r) data.frame(n_sig_snps=nrow(filter(daner_gw, CHR == r$Chr & r$range.left <= BP & BP <= r$range.right))))

region_snp_multiple <- region_snp_counts %>% filter(n_sig_snps > 1)
region_snp_singleton <- region_snp_counts %>% filter(n_sig_snps == 1)

# enumerate regions
region_ranges <- region_snp_multiple %>% mutate(region=row_number())
region_singleton_ranges <- region_snp_singleton %>% mutate(region=row_number()+max(region_ranges$region))

logging(glue("Singleton regions: {nrow(region_snp_singleton)}"))

# clumped SNPs and LD friends
clumped_friends <- bind_rows(plyr::alply(clumped_gw, 1, function(x) tibble(SNP=x$SNP, friend=str_split(x$`LD-friends(0.1).p0.001`, pattern=',')[[1]])))

friends_snp_ld <- str_match(clumped_friends$friend, pattern='(.+)\\((.+)/(.+)\\)')

clumped_friends_ld <- tibble(SNP=clumped_friends$SNP, friend=friends_snp_ld[,2], LD=as.numeric(friends_snp_ld[,3]), kb=as.numeric(friends_snp_ld[,4]))

# Top SNP not in reference panel but LD friend is in COJO results
clumped_friends_cojo <- clumped_friends_ld %>% filter(SNP %in% ref_miss_clump$SNP) %>% filter(friend %in% cojo$SNP)

# daner entries for each COJO SNP
daner_cojo <-
daner %>% 
filter(SNP %in% c(cojo$SNP, ref_miss_cojo_clump$SNP)) %>%
left_join(cojo, by='SNP') 

logging(glue("COJO+Clump SNPs: {nrow(daner_cojo)}"))

# intersection selected SNPs with non-singleton or singleton regions

merge_daner_cojo_regions <- function(daner_cojo, region_ranges) {
    # create genomic range objects for SNPs and regions
    daner_cojo_gr <- with(daner_cojo, GRanges(seqnames=CHR, ranges=IRanges(BP, width=1), SNP=SNP))
    regions_gr <- with(region_ranges, GRanges(seqnames=Chr, ranges=IRanges(start=range.left, end=range.right)))
    
    # find which region each SNP belongs to
    daner_cojo_regions_overlap <- findOverlaps(daner_cojo_gr, regions_gr)
    
    # use overlap object to get row numbers, then slice and bind the two tibbles
    daner_cojo_regions <- 
    bind_cols(dplyr::slice(daner_cojo, daner_cojo_regions_overlap@from),
              dplyr::slice(region_ranges, daner_cojo_regions_overlap@to)) %>%
    # add in a SNP idx if the region cojo file was empty
    mutate(snp_idx=coalesce(snp_idx, 1)) %>%
    # select daner and cojo columns
    select(region, snp_idx, CHR, SNP, BP, A1, A2,
          starts_with('FRQ'), INFO, OR, SE, P, ngt,
          Direction, HetISqt, HetDf, HetPVa, Nca, Nco, Neff_half,
          bJ, bJ_se, pJ, LD_r,
          range.left, range.right,  n_sig_snps) %>%
    arrange(region, snp_idx)
    
    return(daner_cojo_regions)
}

daner_cojo_regions <- merge_daner_cojo_regions(daner_cojo, region_ranges)

daner_cojo_singletons <- merge_daner_cojo_regions(daner_cojo, region_singleton_ranges)

logging(glue("COJO Final SNPs: {nrow(daner_cojo_regions)}"))
logging(glue("COJO Final SNPs p <= 5e-8, pJ <= 5e-8: {nrow(filter(daner_cojo_regions, P <= 5e-8, pJ <= 5e-8))}"))
logging(glue("COJO Final SNPs p > 5e-8, pJ <= 5e-8: {nrow(filter(daner_cojo_regions, P > 5e-8, pJ <= 5e-8))}"))

output_cojo <- snakemake@output$cojo
write_tsv(daner_cojo_regions, output_cojo)

output_singletons <- snakemake@output$singletons
write_tsv(daner_cojo_singletons, output_singletons)







