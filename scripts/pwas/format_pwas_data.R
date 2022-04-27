#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--rosmap", action="store", default=NA, type='character',
              help="Path to ROSMAP data [required]"),
  make_option("--banner", action="store", default=NA, type='character',
              help="Path to ROSMAP data [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

########
# ROSMAP
########

dir.create('resources/data/rosmap_twas/')
unzip(opt$rosmap, exdir='resources/data/rosmap_twas/') 

pos<-fread('resources/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS/train_weights.pos')
pos$N<-376
fwrite(pos, 'resources/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS/train_weights_withN.pos', quote=F, sep=' ', na='NA')

########
# Banner
########

dir.create('resources/data/banner_twas/')
unzip(opt$banner, exdir='resources/data/banner_twas/') 

pos<-fread('resources/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights.pos')
pos$N<-152
fwrite(pos, 'resources/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos', quote=F, sep=' ', na='NA')


