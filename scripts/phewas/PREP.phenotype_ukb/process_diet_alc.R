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


# Load data ---------------------------------------------------------------

fields.all=read.csv(f.,header=T,stringsAsFactors=F)
fields.diet = fields.all %>% .[grep('Diet',.$Path),]

diet=readRDS(f.diet_ac)
diet.online=readRDS(f.diet_online)


# Process -----------------------------------------------------------------

typical_diet_online_label=100020

# Remove fields for pilot data and those that contain texts and types
remove_category <- function(x,kw,col.n,exact=F){
  if (exact==T){
    new.x=x[!grepl(paste0('^',kw,'$'),x[,col.n]),]
  }else{
    new.x=x[!grepl(kw,x[,col.n]),]
  }
  return(new.x)
}

fields.diet = fields.diet %>%
  remove_category(.,kw='(pilot)',col.n='Field') %>%
  remove_category(.,kw='Reason',col.n='Field') %>%
  (function(x) x[,!colnames(x) %in% c(100020,20085,20086)]) %>%
  remove_category(.,kw='type',col.n='Field') %>%
  remove_category(.,kw='Type',col.n='Field')

