library(data.table)
library(dplyr)
library(knitr)
library(TwoSampleMR)
library(MRPRESSO)
library(ggplot2)
library(optparse)
library(ggpubr)
library(readr)
options(bitmapType='cairo') # Specific to Eddie - to enable writing figures

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

option_list <- list(
  make_option('--MRexposure', type='character', help="Intermediate files of exposure data for MR", action='store'),
  make_option('--MRoutcome', type='character', help='Intermediate files of outcome data for MR', action='store'),
  make_option('--outcomeName', default='MDD', type='character', help='Name of outcome phenotype', action='store'),
  make_option('--outputfig', type='character', help='Folder for QC figures', action='store'),
  make_option('--outputtable', type='character', help='Summary file for MR analysis', action='store'),
  make_option('--saveHarmonisedData', type='character', help='Save harmonised data for analysis', action='store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

f.exposure=opt$MRexposure
f.outcome=opt$MRoutcome
outcome.name=opt$outcomeName
ofig.path=opt$outputfig
otable.path=opt$outputtable
saveDat=opt$saveHarmonisedData

# Basic settings ----------------------------------------------------------
dir.create(file.path(ofig.path), showWarnings = FALSE)
#dir.create(file.path(otable.path), showWarnings = FALSE)


# MR analysis -------------------------------------------------------------

  ## exposure dat
  exposure_dat=read_tsv(f.exposure)
  ## outcome dat
  outcome_dat=read_tsv(file=f.outcome)
  
  # Harmonise the exposure and outcome data
  dat <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  res <- mv_multiple(dat,plots = T)
  
  res$result %>% as.data.frame %>%
  write_tsv(.,file=otable.path,
            col_names = T)
  
  names(res$plots) = res$result$exposure
  names(res$plots) %>% as.list %>% 
    lapply(.,FUN=function(x) ggsave(res$plots[[x]],file=paste0(ofig.path,'/',x,'_plots.png'), 
                    width=10, height=15,dpi=300))

rm(list=ls())