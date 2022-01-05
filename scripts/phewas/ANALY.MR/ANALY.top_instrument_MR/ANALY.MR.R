library(dplyr)
library(readr)
library(knitr)
library(TwoSampleMR)
library(MRPRESSO)
library(ggplot2)
library(optparse)
library(ggpubr)
library(pbapply)
library(harrypotter)
options(bitmapType='cairo') # Specific to Eddie - to enable writing figures


# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1

parse <- OptionParser()

option_list <- list(
  make_option('--ListMR', type='character', help="A list of MR analysis to run", action='store'),
  make_option('--out', type='character', help="Path for outputs", action='store'),
  make_option('--saveHarmonisedData', default=T,type='logical', help='Save harmonised data for analysis?', action='store'),
  make_option('--NoSplit', default=T,type='logical', help='Save harmonised data for analysis?', action='store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

d.mr_input = opt$ListMR
d.output = opt$out
do.saveHarmoDat = opt$saveHarmonisedData
do.no_split = opt$NoSplit

# d.mr_input = 'data/MR/ls.mr_analysis.rds'
# d.output = 'results/phewas/MR/ukb'
# do.saveHarmoDat = T
# do.no_split = T

# Basic settings ----------------------------------------------------------

dir.create(file.path(d.output), showWarnings = FALSE)

# Load inputs -------------------------------------------------------------

inputs.mr = readRDS(d.mr_input)

# Define function for analysis --------------------------------------------

batch_mr <- function(x.inputs){
  
  f.exposure = x.inputs$file_exposure
  f.outcome = x.inputs$file_outcome
  prefix.out = x.inputs$output_prefix
  info.input = x.inputs %>% select(-file_exposure,-file_outcome,-output_prefix)
  
  # Load data
  dat.exposure = read_tsv(f.exposure)
  dat.outcome = read_tsv(f.outcome)
  
  dat.mr = harmonise_data(dat.exposure,dat.outcome) %>% 
    .[order(.$pval.exposure,decreasing = F),]
  if(nrow(dat.mr)>=121&(dat.mr$exposure=='MDD')){dat.mr = dat.mr %>% .[1:121,]}else{}
  if (do.saveHarmoDat==T){      
    write_tsv(dat.mr,paste0(d.output,'/',prefix.out,"_mrDat.csv"))
  }
  
  # MR analysis
  MR=mr(dat = dat.mr,
        method_list =c("mr_ivw",
                       "mr_egger_regression",
                       "mr_egger_regression_bootstrap",
                       "mr_weighted_median",
                       "mr_sign"))
  het = mr_heterogeneity(dat.mr) %>% 
    filter(method=='Inverse variance weighted') %>% 
    select(Q,Q_df,Q_pval)
  plt= mr_pleiotropy_test(dat.mr) %>% 
    select(egger_intercept,egger_se=se,egger_pval=pval)
  res_single = mr_singlesnp(dat.mr)
  res_loo = mr_leaveoneout(dat.mr)

  
  # Summarise results
  mr_res=data.frame(MR,plt,het,info.input, stringsAsFactors = F)
  write_tsv(mr_res,paste0(d.output,'/mr_res.',prefix.out,'.tsv'))
  
# Create plots
  fig.scatter = mr_scatter_plot(MR, dat.mr)[[1]]+scale_colour_hp_d(option = "ronweasley2", name = "MR method")
  fig.forest = mr_forest_plot(res_single)[[1]]+scale_colour_hp_d(option = "ronweasley2", name = "MR method")
  fig.loo = mr_leaveoneout_plot(res_loo)[[1]]
  fig.funnel = mr_funnel_plot(res_single)[[1]]
  
  # Save plots per test
  fig.total = ggarrange(fig.scatter,fig.funnel,
                        fig.forest,fig.loo,
                        ncol = 2,nrow=2,
                        labels = c('Scatter plot','Funnel plot',
                                   'Forest plot','Leave-one-out plot'),
                        widths = c(1,1),heights = c(0.5,1.5),align = 'v')
  ggsave(fig.total, file=paste0(d.output,'/',prefix.out,"_plot.png"), 
         width=10, height=15,dpi=300)
  
  # Save scatter plots for further processing
  fig.scatter = fig.scatter
  save(fig.scatter,file=paste0(d.output,'/',prefix.out,"_scatterplot.RData"))
  
  gc()
}


# Run MR ------------------------------------------------------------------

inputs.mr %>% split(.,seq(nrow(.))) %>% 
  pblapply(batch_mr)


# Put all results together ------------------------------------------------

if (do.no_split==T){
  ls.f.res = list.files(path=d.output,pattern = 'mr_res.',full.names = T)
  all.res = ls.f.res %>% as.list %>% 
    pblapply(read_tsv) %>% 
    bind_rows
  write_tsv(all.res,paste0(d.output,'/ALL_mr_res.tsv'))
  system(paste0('rm ',d.output,'/mr_res*'))
}
