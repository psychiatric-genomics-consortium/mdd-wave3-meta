library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dic = args[1]
f.diet_ac = args[2]                 # file location for cognition results (assessment centre)
f.diet_online = args[3]             # file location for cognition results (online)
f.cogvar_QCed = args[4]            # file location for QCed variables (non-outcome vars removed, check Extract field list section)
f.output_data = args[5]
f.output_dictionary = args[6]
