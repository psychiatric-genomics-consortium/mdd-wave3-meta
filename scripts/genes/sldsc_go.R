# Merge results from S-LDSC analysis of GO pathways

library(dplyr)
library(readr)
library(stringr)

save.image()

# get list of S-LDSC results files and name them with the geneset
sldsc_results <- snakemake@input$sldsc
names(sldsc_results) <- snakemake@params$genesets

# read in results files
sldsc <- lapply(sldsc_results, read_tsv)

# extract enrichment for genesets, which are the final row of
# each results file with category "L2_1"
# Calculate FDR-adjusted P
geneset_sldsc <- bind_rows(lapply(sldsc, function(result) filter(result, Category == 'L2_1')), .id='GeneSet') %>%
mutate(Enrichment_p_fdr=p.adjust(Enrichment_p, 'fdr'))

write_tsv(geneset_sldsc, snakemake@output[[1]])