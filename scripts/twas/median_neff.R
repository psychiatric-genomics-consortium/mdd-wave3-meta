#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--munged", action="store", default=NA, type='character',
              help="Name of munged sumstats file [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Name of output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
ss<-fread(opt$munged)
med_N<-median(ss$N, na.rm=T)
write.table(med_N, opt$out, col.names=F, row.names=F, quote=F)

