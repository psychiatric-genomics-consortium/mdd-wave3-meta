
### Strings aren't factors

options(stringsAsFactors=F)

### Initialise args

args <- commandArgs(trailingOnly=TRUE)

### Load packages

library(data.table)

## Define input and output

input <- args[1]
output <- args[2]

## Load COJO data

COJO <- fread(input)

## Initialise start

start <- 0

## Initialise output dataframe

Outputdf <- data.frame(Loci = character(length=dim(COJO)[1]))

for(i in 1:dim(COJO)[1]){
    ## Set start - if the COJO SNP is within the first (second) Mb, set to 1 (1000001), otherwise set start so that the COJO SNP is approximately central
    if(floor(COJO$BP[i] / 1000000) == 0){
        start <- 1 
    } else if(floor(COJO$BP[i] / 1000000) == 1){
        start <- 1000001
    } else
    start <- 1 + ((floor(COJO$BP[i] / 1000000)-1) * 1000000)
    Outputdf$Loci[i]<-as.character(paste(COJO$CHR[i], ".", start, "_", start + 3000000, sep=""))
}

## Write output

fwrite(file=output, Outputdf, col.names=F)
