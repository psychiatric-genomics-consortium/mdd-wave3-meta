library(dplyr)
library(data.table)
library(readr)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dictionary = args[1]
f.experiment_procedure = args[2]
f.covered_by_sumscores = args[3]
f.output = args[4]



# Load data ---------------------------------------------------------------
dat.dic=read.csv(f.dictionary,header=T,stringsAsFactors=F)

