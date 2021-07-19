# Merge COJO regions together and compare to clumped results,
# regions plots to create list of independent SNPs

library(dplyr)
library(readr)
library(stringr)
library(glue)

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
    cojo_stats <- read_tsv(cojo_file, col_types=cols('SNP'=col_character(), 'refA'=col_character()))
    # parse region coordinates from filename
    regions_matched = str_match(basename(cojo_file), "([[:digit:]]+):([[:digit:]]+)-([[:digit:]]+)")
    cojo_region <- cojo_stats %>%
    mutate(region.left=as.numeric(regions_matched[,3]),
           region.right=as.numeric(regions_matched[,4]),
           ridx=row_number())
    return(cojo_region)
}


cojo_list <- lapply(regions_jma_cojo, read_cojo)
cojo <- bind_rows(cojo_list) 

# region coordinates
regions <- cojo %>% group_by(Chr, region.left, region.right) %>% summarise(p=min(p)) %>% arrange(p) %>% ungroup() 

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

# count number of significant SNPs in each region
region_snp_counts <- plyr::adply(regions, 1, function(r) data.frame(n_sig_snps=nrow(filter(daner_gw, CHR == r$Chr & r$region.left <= BP & BP <= r$region.right)))) %>% arrange(Chr, region.left)

region_snp_multiple <- region_snp_counts %>% filter(n_sig_snps > 1)
region_snp_singleton <- region_snp_counts %>% filter(n_sig_snps == 1)

logging(glue("Singleton regions: {nrow(region_snp_singleton)}"))

# enumerate regions
region_numbers <- region_snp_multiple %>% arrange(p) %>% mutate(region=row_number())

# clumped SNPs and LD friends
clumped_friends <- bind_rows(plyr::alply(clumped_gw, 1, function(x) tibble(SNP=x$SNP, friend=str_split(x$`LD-friends(0.1).p0.001`, pattern=',')[[1]])))

friends_snp_ld <- str_match(clumped_friends$friend, pattern='(.+)\\((.+)/(.+)\\)')

clumped_friends_ld <- tibble(SNP=clumped_friends$SNP, friend=friends_snp_ld[,2], LD=as.numeric(friends_snp_ld[,3]), kb=as.numeric(friends_snp_ld[,4]))

# Top SNP not in reference panel but LD friend is in COJO results
clumped_friends_cojo <- clumped_friends_ld %>% filter(SNP %in% ref_miss_clump$SNP) %>% filter(friend %in% cojo$SNP)

# daner entries for each COJO SNP, removing singletons
daner_cojo <-
daner %>% 
inner_join(cojo, by='SNP') %>%
inner_join(select(region_numbers, 'Chr', 'region', 'region.left'), by=c('CHR'='Chr', 'region.left')) %>%
mutate(locus=rank(P)) %>%
select(region, ridx, locus, CHR, SNP, A1, A2, starts_with('FRQ'), INFO, OR, SE, P, ngt,
      Direction, HetISqt, HetDf, HetPVa, Nca, Nco, Neff_half, bJ, bJ_se, pJ, LD_r, region.left, region.right) %>%
arrange(region, ridx)

logging(glue("COJO Final SNPs: {nrow(daner_cojo)}"))
logging(glue("COJO Final SNPs p <= 5e-8, pJ <= 5e-8: {nrow(filter(daner_cojo, P <= 5e-8, pJ <= 5e-8))}"))
logging(glue("COJO Final SNPs p > 5e-8, pJ <= 5e-8: {nrow(filter(daner_cojo, P > 5e-8, pJ <= 5e-8))}"))

output_cojo <- snakemake@output[[1]]

write_tsv(daner_cojo, output_cojo)







