library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(optparse)

source('scripts/phewas/util/reg_phewasStyle.R')

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

args = commandArgs(trailingOnly=TRUE)

f.dat_dir = args[1]             
f.PRS = args[2]                 
f.model = args[3]
target.category = args[4]
#f.output_res = args[5]

# f.dat_dir = 'data/dat.imaging_chunk.rds'
# f.PRS = 'data/PRS_all.rds'
# f.model = 'results/phewas/models.rds'
# f.output_res = 'results/phewas/phewas_out.rds'



# Load phenotyp -----------------------------------------------------------

g.path = f.dat_dir %>% strsplit(.,'/') %>% unlist %>% 
  .[1:(length(.)-1)] %>% 
  paste0(.,collapse = '/')

fs = list.files(path=g.path,pattern = 'dat.', full.names = T)%>% 
  (function(x) x[grep('chunk.rds',x)])
phewas.dat = fs %>% as.list %>% lapply(readRDS) %>% 
  Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="f.eid",all=T), .) %>% 
  right_join(.,readRDS(f.PRS),by='f.eid')

ls.models = readRDS(f.model)

if(target.category!='all'){
  ls.models = ls.models %>% filter(category==target.category)
}else if(target.category=='all'){
  ls.models = ls.models
}

# Remove phenotypes that do not have enough variance for multiple assessment centres
ls.models = ls.models 
#%>% .[!.$dependent %in% c('f.6153','f.6177','f.2375','f.2385','f.2395','f.2405','f.2674','f.2684','f.2694','f.2704','f.2714','f.2724','f.2734','f.2744','f.2754','f.2764','f.2774','f.2784','f.2794','f.2804','f.2814','f.2824','f.2834','f.3140','f.3536','f.3581','f.3591','f.3700','f.3710','f.3720','f.3829','f.3839','f.3849','f.3872','f.3882','f.4041','f.3546'),]


result=reg_phewasStyle(ls.models,dat_short=phewas.dat)

target.category = target.category %>% 
  gsub(' ','_',.)
target.category
saveRDS(result,paste0('results/phewas/phewas_out_',target.category,'.rds'))


