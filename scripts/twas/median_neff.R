#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--daner", action="store", default=NA, type='character',
              help="Name of daner file [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Name of output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

if(substr(opt$daner,(nchar(opt$daner)+1)-3,nchar(opt$daner)) == '.gz'){
  ss<-fread(cmd=paste0('zcat ',opt$daner), nrow=1)
} else {	
  ss<-fread(opt$daner, nrow=1)
}

Neff_col<-which(names(ss) == 'Neff_half')

if(substr(opt$daner,(nchar(opt$daner)+1)-3,nchar(opt$daner)) == '.gz'){
  ss<-fread(cmd=paste0('zcat ',opt$daner,' | cut -f 19'))
} else {	
  ss<-fread(cmd=paste0('cut -f 19 ',opt$daner))
}

Neff<-as.numeric(ss[['Neff_half']])
med_Neff<-median(Neff, na.rm=T)

write.table(med_Neff, opt$out, col.names=F, row.names=F, quote=F)

