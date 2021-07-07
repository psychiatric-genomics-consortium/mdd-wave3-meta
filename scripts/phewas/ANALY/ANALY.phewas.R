# Basic setups ------------------------------------------------------------

library('dplyr')
library('pbapply')
library('nlme')

setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/')

phewas.dat=readRDS('data/dat.phewas_imaging.rds')

# Define functions --------------------------------------------------------
source('FUNs/reg_phewasStyle.R')

targetdata=targetdata
