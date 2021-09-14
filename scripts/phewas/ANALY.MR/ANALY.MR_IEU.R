library(dplyr)
library(readr)
library(knitr)
library(TwoSampleMR)
library(MRPRESSO)
library(ggplot2)
library(optparse)
library(ggpubr)
library(pbapply)
options(bitmapType='cairo') # Specific to Eddie - to enable writing figures


# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1

parse <- OptionParser()

option_list <- list(
  make_option('--ListPheno', type='character', help="A list of IEU IDs to analyse", action='store'),
  make_option('--MDDsumstats', type='character', help="MDD full sumstats", action='store'),
  make_option('--MDDinstrument', type='character', help="MDD instruments", action='store'),
  make_option('--out', type='character', help="Path for outputs", action='store'),
  make_option('--saveHarmonisedData', default=T,type='logical', help='Save harmonised data for analysis?', action='store'),
  make_option('--NoSplit', default=T,type='logical', help='Save harmonised data for analysis?', action='store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

f.mr_input = opt$ListMR
f.mdd_sumstats = opt$MDDsumstats
f.mdd_instruments = opt$MDDinstrument
d.output = opt$out
do.saveHarmoDat = opt$SaveHarmonisedData
do.no_split = opt$NoSplit

# f.mr_input = 'docs/tables/ldsc_open_mr_candidates.txt'
# f.mdd_sumstats = 'results/distribution/daner_pgc_mdd_noUKBB_eur_hg19_v3.49.24.05.rp.gz'
# f.mdd_instruments = 'data/MR/MR_InterFiles/mdd_exposure_outcome/MDD.exposure_dat'
# d.output = 'results/phewas/MR/mr_base'
# do.saveHarmoDat = T
# do.no_split = T

# Basic settings ----------------------------------------------------------

dir.create(file.path(d.output), showWarnings = FALSE)

# Load inputs -------------------------------------------------------------

inputs.mr = read_tsv(f.mr_input) %>% 
  .[!.$id %in% c('ebi-a-GCST005047'),]

mdd.sumstats = read_tsv(f.mdd_sumstats) %>% 
  mutate(N=2*Neff_half,Effect=log(OR)) %>% 
  select(SNP,A1,A2,Effect,SE,P,Direction,af=FRQ_A_235237,CHR,BP,N) 

mdd.exposure = read_tsv(f.mdd_instruments)

# Define function for analysis --------------------------------------------

batch_mr <- function(x.input,mdd_type){
  
  id.pheno = x.input$id
  if(mdd_type=='exposure'){
    prefix.out=paste0('MDD_',id.pheno) 
  }else if(mdd_type=='outcome'){
    prefix.out=paste0(id.pheno,'_MDD')
  }
  
  info.input = x.input %>% 
    select(trait,rg,rg_p=p,rg_qval=qvalue,gcov_int,subcategory) 
  
  # Load data
  if(mdd_type=='exposure'){
    dat.exposure = mdd.exposure
    dat.outcome = extract_outcome_data(snps = dat.exposure$SNP, outcomes = id.pheno)
  }else if(mdd_type=='outcome'){
    dat.exposure = extract_instruments(outcomes=id.pheno)
    if(!is.null(nrow(dat.exposure))){
      match.outcome = mdd.sumstats %>% .[.$SNP %in% dat.exposure$SNP,] %>%
        dplyr::select(SNP,effect_allele=A1,other_allele=A2,
                      eaf=af,beta=Effect,se=SE,pval=P,samplesize=N) %>% 
        mutate(Phenotype='MDD')
      if(!is.null(nrow(match.outcome))){
        if(nrow(match.outcome)>10){
        dat.outcome = format_data(match.outcome, type="outcome")
        }
      }
    }
  }
  
  if(exists('dat.outcome')){
    dat.mr = harmonise_data(dat.exposure,dat.outcome)
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
    fig.scatter = mr_scatter_plot(MR, dat.mr)
    fig.forest = mr_forest_plot(res_single)
    fig.loo = mr_leaveoneout_plot(res_loo)
    fig.funnel = mr_funnel_plot(res_single)
    
    # Save plots per test
    fig.total = ggarrange(fig.scatter[[1]],fig.funnel[[1]],
                          fig.forest[[1]],fig.loo[[1]],
                          ncol = 2,nrow=2,
                          labels = c('Scatter plot','Funnel plot',
                                     'Forest plot','Leave-one-out plot'),
                          widths = c(1,1),heights = c(1,1.5),align = 'v')
    ggsave(fig.total, file=paste0(d.output,'/',prefix.out,"_plot.png"), 
           width=10, height=15,dpi=300)
    
    # Save scatter plots for further processing
    fig.scatter = fig.scatter[[1]]
    save(fig.scatter,file=paste0(d.output,'/',prefix.out,"_scatterplot.RData"))
  }
  gc()
}


# Run MR ------------------------------------------------------------------

inputs.mr %>% split(.,seq(nrow(.))) %>% 
  pblapply(batch_mr,mdd_type='exposure')
inputs.mr %>% split(.,seq(nrow(.))) %>% 
  pblapply(batch_mr,mdd_type='outcome')

# Put all results together ------------------------------------------------

if (do.no_split==T){
  ls.f.res = list.files(path=d.output,pattern = 'mr_res.',full.names = T)
  all.res = ls.f.res %>% as.list %>% 
    pblapply(read_tsv) %>% 
    bind_rows
  write_tsv(all.res,paste0(d.output,'/ALL_mr_res.tsv'))
  system(paste0('rm ',d.output,'/mr_res*'))
}
