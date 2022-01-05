library(dplyr)
library(data.table)
library(readr)
library(TwoSampleMR)
library(pbapply)
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1

parse <- OptionParser()

option_list <- list(
   make_option('--gwas', type='character', help="GWAS sumstats", action='store'),
   make_option('--ExposureProcess', type='character', default = T, help="Process as exposure", action='store'),
   make_option('--OutcomeProcess', type='character', default = T, help="Process as outcome", action='store'),
   make_option('--ExposureDat', type='character', default = "",help="Location of exposure data (OutcomeProcess==T)", action='store'),
   make_option('--fNeale', type='character', default = "data/MR/MR_sumstats/Neale/variants.tsv.gz",
               help="Location of info file for Neale lab", action='store'),
   make_option('--fBIG40', type='character', default = "data/MR/MR_sumstats/BIG40/variants.txt.gz",
               help="Location of info file for Neale lab", action='store'),
   make_option('--phenoSummary', type='character', default = "data/MR/pheno_gwas_n.rds",
               help="Location of info file for Neale lab", action='store'),
   make_option('--interexp', type='character', 
               help='Folder for intermediate files', action='store')
)

args = commandArgs(trailingOnly=TRUE)

opt <- parse_args(OptionParser(option_list=option_list), args=args)
gwas.f_and_label = opt$gwas
f.exposure = opt$ExposureDat
d.output=opt$interexp

do.exposure = opt$ExposureProcess %>% as.logical
do.outcome = opt$OutcomeProcess %>% as.logical
f.Neale = opt$fNeal
f.big40 = opt$fBIG40
f.nsummary = opt$phenoSummary

# gwas.f_and_label= 'data/MR/pheno_gwas_forMR.rds'
# f.exposure='data/MR/MR_InterFiles/exposure_stats/MDD.exposure_dat'
# d.output='data/MR/MR_InterFiles'
# do.exposure=T
# do.outcome=T
# f.Neale='data/MR/MR_sumstats/Neale/variants.tsv.gz'
# f.big40='data/MR/MR_sumstats/BIG40/variants.txt.gz'
# f.nsummary='data/MR/pheno_gwas_n.rds'


# Basic settings ----------------------------------------------------------

dir.create(file.path(d.output), showWarnings = FALSE)

# inputs
inputs.gwas = readRDS(gwas.f_and_label)

### Load ref data
# Neale gwas
ref.neale = fread(f.Neale) %>% 
   filter(info>0.8,minor_AF>0.0005) %>% 
   select(variant,A1=ref,A2=alt,rsid,info,freq=AF,chr,pos,n=n_called)
# BIG40 IDP gwas
ref.big40 = fread(f.big40) 
# Phenotype N
ref.n = readRDS(f.nsummary)

if(do.outcome==T){
   external.exposure = read_tsv(f.exposure) %>% .$SNP
}


# Define process ----------------------------------------------------------

process_gwas <- function(tmp.input,tmp_process.exposure,tmp_process.outcome,ls.snp=NA){
   tmp.fname = tmp.input$file_loc
   target.colname = tmp.input$field_tag
   tmp.type = tmp.input$loading_type
   
   
   target.pheno = target.colname %>% strsplit(.,'.',fixed = T) %>% unlist %>% head(n=2) %>%
      tail(n=1)
   
   tmp.n = ref.n$n[ref.n$phenotype==target.pheno] %>% .[1]
   
   
   ################ Create a list that do not need to process #####################
   ls.processed = list.files(d.output,pattern = '_dat$',full.names = T,recursive = T) %>% 
      basename %>% 
      gsub('.exposure_dat','',.,fixed = T) %>% 
      gsub('.outcome_dat','',.,fixed = T)
   
   ls.running = read_tsv('data/MR/running_prepare_expo_outc',col_names = F) %>% as.data.frame %>% 
         .$X1 %>% as.character %>% as.vector
   
   if(file.exists(paste0(d.output,'/Instrument_info.tsv'))){
      ls.no_instrument = read_tsv(paste0(d.output,'/Instrument_info.tsv'),col_names = F) %>% 
         as.data.frame %>%  
         .$X1 
   }else{ls.no_instrument=NA}
   
   ls.avoid = c(ls.processed,ls.running,ls.no_instrument)
   
   ls.current = target.colname
   ################################################################################ 
   
   
   if (sum(ls.current %in% ls.avoid)==0){
      cat('\nLoading gwas sumstats\n')
      system(paste0('echo ',ls.current,' >> data/MR/running_prepare_expo_outc'))
   
      if (tmp.type=='neale'){
         gwas_sumstats = read_tsv(tmp.fname) %>% 
            left_join(.,ref.neale,by='variant') %>% 
            filter(chr!='X') %>%
            mutate(chr=as.numeric(chr)) %>% 
            select(SNP=rsid,Allele1=A2,Allele2=A1,
                   Effect=beta,se,p=pval,Freq1=freq,CHR=chr,BP=pos) %>% 
            mutate(N=tmp.n) 
         
      }else if(tmp.type=='big40'){
         gwas_sumstats = read_delim(tmp.fname,delim=" ")%>% 
            mutate(pval=10^(-`pval(-log10)`),n=tmp.n) %>%
            left_join(.,ref.big40,by=c('chr'='chr','pos'='pos','rsid'='rsid','a1'='a1','a2'='a2')) %>% 
            filter(info>0.8,chr!='0X') %>% 
            mutate(chr=as.numeric(chr)) %>% 
            select(SNP=rsid,Allele1=a2,Allele2=a1,
                   Effect=beta,se,p=pval,Freq1=af,CHR=chr,BP=pos) %>% 
            mutate(N=tmp.n) %>% 
            as.data.frame
      }else if(tmp.type=='mtag'){
         gwas_sumstats = read_tsv(tmp.fname)%>% 
            left_join(.,ref.big40,by=c('BP'='pos','SNP'='rsid','A1'='a1','A2'='a2')) %>% 
            select(SNP,Allele1=A1,Allele2=A2,
                   Effect=mtag_beta,se=mtag_se,p=mtag_pval,Freq1=meta_freq,CHR,BP) %>%
            mutate(N=tmp.n) %>% 
            as.data.frame
      }else if(tmp.type=='local'){
         gwas_sumstats = read_tsv(tmp.fname) %>% 
            as.data.frame
      }else{
         print('No loading type')
      }
      cat('Sumstats loaded\n')
   }
   
   if (tmp_process.exposure==T & sum(ls.current %in% ls.avoid)==0){
      cat('Processing as exposure\n')
      ### find top SNPs
      SigniSNP=gwas_sumstats %>% 
         mutate(NewAllele1=ifelse(Effect>=0,Allele1,Allele2),
                NewAllele2=ifelse(Effect>=0,Allele2,Allele1),
                NewFreq1=ifelse(Effect>=0,Freq1,1-Freq1),
                NewEffect=ifelse(Effect>=0,Effect,-1*Effect)) %>% 
         as.data.frame %>% 
         select(SNP,effect_allele=NewAllele1,other_allele=NewAllele2,
                eaf=NewFreq1,beta=NewEffect,se,pval=p,samplesize=N) %>% 
         mutate(units='unit',Phenotype=target.colname) %>% 
         filter(pval<5e-8)
      
      ### clump exposure dat
      if(nrow(SigniSNP)==0){
         info.summary = data.frame(phenotype=target.colname,
                                   n_5e_08=0,
                                   n_1e_07=0,
                                   n_1e_06=0,
                                   ninstruments=0,
                                   dat_loc=NA)
      }else{
         exposure_dat = format_data(SigniSNP, type="exposure")
         exposure_dat <- clump_data(exposure_dat,clump_kb=500,clump_r2=0.001)
         
         info.summary = data.frame(phenotype=target.colname,
                                   n_5e_08=sum(gwas_sumstats$p<5e-8),
                                   n_1e_07=sum(gwas_sumstats$p<1e-7),
                                   n_1e_06=sum(gwas_sumstats$p<1e-6),
                                   ninstruments=nrow(exposure_dat),
                                   dat_loc=NA)

         if(nrow(exposure_dat)>=10){
            write_tsv(exposure_dat,paste0(d.output,'/',target.colname,".exposure_dat"))
            info.summary$dat_loc=paste0(d.output,'/',target.colname,".exposure_dat")
         }
      }
      write_tsv(info.summary,paste0(d.output,'/Instrument_info.tsv'),col_names = F,append = T)
   }else{
      cat('Exposure processing skipped\n')
   }
   
   if (tmp_process.outcome==T & tmp_process.exposure==F & sum(ls.current %in% ls.avoid)==0){
      cat('Processing as outcome\n')
      outcome_dat=gwas_sumstats %>% 
         as.data.frame %>% 
         select(SNP,effect_allele=Allele1,other_allele=Allele2,
                eaf=Freq1,beta=Effect,se,pval=p,samplesize=N) %>% 
         mutate(units='unit',Phenotype=target.colname) %>% 
         .[.$SNP %in% ls.snp,] 
      outcome_dat = format_data(outcome_dat, type="outcome")
      write_tsv(outcome_dat,paste0(d.output,'/',target.colname,".outcome_dat"))
   }else if(tmp_process.outcome==T & tmp_process.exposure==T & exists('exposure_dat') & sum(ls.current %in% ls.avoid)==0){
      if(nrow(exposure_dat)>=10){
         cat('Processing as outcome\n')
         outcome_dat=gwas_sumstats %>% 
            as.data.frame %>% 
            select(SNP,effect_allele=Allele1,other_allele=Allele2,
                   eaf=Freq1,beta=Effect,se,pval=p,samplesize=N) %>% 
            mutate(units='unit',Phenotype=target.colname) %>% 
            .[.$SNP %in% ls.snp,] 
         outcome_dat = format_data(outcome_dat, type="outcome")
         write_tsv(outcome_dat,paste0(d.output,'/',target.colname,".outcome_dat"))
      }
   }else{
      cat('Outcome processing skipped\n')
   }
   gc()
}


# Process exposure data ---------------------------------------------------

ls.rm.1 = list.files(d.output,pattern = '_dat$',full.names = T,recursive = T) %>% 
   basename %>% 
   gsub('.exposure_dat','',.,fixed = T) %>% 
   gsub('.outcome_dat','',.,fixed = T)
if(file.exists(paste0(d.output,'/Instrument_info.tsv'))){
   ls.rm.2 = read_tsv(paste0(d.output,'/Instrument_info.tsv'),col_names = F) %>% 
      as.data.frame %>%  
      .$X1 
}else{ls.rm.2=NA}
ls.rm = c(ls.rm.1,ls.rm.2)

system('touch data/MR/running_prepare_expo_outc')
inputs.gwas %>% .[!.$field_tag %in% ls.rm,]%>% split(.,seq(nrow(.))) %>%
   pblapply(.,process_gwas,
            tmp_process.exposure=do.exposure,
            tmp_process.outcome=do.outcome,ls.snp=external.exposure)

rm(list=ls())