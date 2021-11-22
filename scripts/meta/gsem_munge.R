library(GenomicSEM)

munge(files=snakemake@input$sumstats,
      hm3=snakemake@input$hm3,
      trait.names=snakemake@params$prefix,
      info.filter=0.9,
      maf.filter=0.01)