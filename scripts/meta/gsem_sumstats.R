# Create sumstats for user GWAS of phenotype-specific meta-analyses

library(GenomicSEM)
library(readr)


mdd_sumstats <- sumstats(files=snakemake@input$sumstats,
                         ref=snakemake@input$ref,
                         trait.names=snakemake@params$cohorts,
                         se.logit=rep(T, times=length(snakemake@params$cohorts)),
                         info.filter=0.6,
                         maf.filter=0.01,
                         parallel=TRUE,
                         cores=8)
                         
write_tsv(mdd_sumstats, snakemake@output[[1]])