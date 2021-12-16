
### Initialise args

args <- commandArgs(trailingOnly=TRUE)

### Load packages

library(data.table)

## Define input and output

input <- args[1]
output <- gsub(gsub(input, pattern=".gz", replacement=""), pattern="finemapping/locus_results", replacement="finemapping/locus_credible_causal")

## Load finemap data

Finemap <- fread(cmd=paste("gunzip -c ", input, sep=""))

## Set up counters

PIP_Tot <- 0
rownum <- 0

## ID credible causal set

while(PIP_Tot < 0.95){
  rownum <- rownum + 1
  PIP_Tot <- PIP_Tot + Finemap$PIP[rownum]
}

## Write output

fwrite(Finemap[1:rownum], output, sep="\t")