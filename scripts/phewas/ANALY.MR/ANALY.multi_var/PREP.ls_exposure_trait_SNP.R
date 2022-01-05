library(dplyr)
library(tidyverse)
library(data.table)
library(readr)
library(pbapply)
library(ggplot2)
library(knitr)
library(kableExtra)
library(here)
library(ggpubr)

source(here::here('scripts/phewas/util/PheWAS_style_p_plot.R'))
source(here::here('scripts/phewas/util/beta_p_plots_sig.R'))


# Load single trait MR results --------------------------------------------

f.category = 'data/phewas_categories.tsv'
ref.category = fread(here::here(f.category),header=F,sep='\n') %>% 
  as.data.frame %>% .$V1


mr.ukb=read_tsv(here::here('results/phewas/MR/ukb/ALL_mr_res.tsv')) %>% 
  .[grepl('^Inverse variance weighted$|^MR Egger$|Weighted median',.$method),] %>% 
  mutate(label_exposure=label_exposure %>% gsub(' (left)','',.,fixed = T) %>% 
           gsub( '(left hemisphere)','',.,fixed = T) %>%
           gsub( '(right hemisphere)','',.,fixed = T) %>%
           gsub(' (whole brain)','',.,fixed = T)) %>% 
  mutate(label_outcome=label_outcome %>% gsub(' (left)','',.,fixed = T) %>% 
           gsub( '(left hemisphere)','',.,fixed = T) %>%
           gsub( '(right hemisphere)','',.,fixed = T) %>%
           gsub(' (whole brain)','',.,fixed = T))

ref.two_level.category = data.frame(lv2=ref.category,stringsAsFactors = F) %>%
  mutate(lv1=ifelse(grepl('Mental health',lv2),'Mental health',
                    ifelse(grepl('Early life|Socio|Lifestyle|Diet|Physical activity',lv2),'Environment',
                           ifelse(grepl('Physical health|Physical|Blood|Body MRI',lv2),'Self-reported and assessed physical condition',
                                  'Brain and cognition'))))

mr.ukb = mr.ukb %>% 
  left_join(.,ref.two_level.category,by=c('category'='lv2')) 



# Find traits to enter multivariable MR -----------------------------------

### Find traits using three methods

find_sig <- function(x.res,x.pheno,target.column){
  x.res = x.res %>% 
    filter(get(target.column)==x.pheno)
  sig.output = ifelse(x.res$p.fdr[x.res$method=='Inverse variance weighted']<0.05,1,0)+
    ifelse(x.res$egger_p.fdr[x.res$method=='MR Egger']>0.05,1,
           ifelse(x.res$p.fdr[x.res$method=='MR Egger']<0.05,1,0))+
    ifelse(x.res$Q_p.fdr[x.res$method=='Weighted median']>0.05,1,
           ifelse(x.res$p.fdr[x.res$method=='Weighted median']<0.05,1,0))+
    ifelse(abs(sum(sign(x.res$b)))==3,1,0)
  output = sig.output==4
  output = data.frame(pheno=unique(x.res[,target.column]),sig.output,is.valid=output)
  return(output)
}



mr.ukb.pheno_to_mdd = mr.ukb %>% filter(outcome=='MDD',nsnp>=30) %>% 
  group_by(method) %>% 
  mutate(p.fdr = p.adjust(pval,method='fdr'),
         Q_p.fdr = p.adjust(Q_pval,method='fdr'),
         egger_p.fdr = p.adjust(egger_pval,method='fdr'))

ls.sig.pheno_to_mdd = mr.ukb.pheno_to_mdd$exposure %>% unique %>% as.list %>% 
  pblapply(.,find_sig,x.res=mr.ukb.pheno_to_mdd,target.column='exposure') %>% 
  bind_rows %>% 
  filter(is.valid==T)

mr.sig.ukb.pheno_to_mdd = mr.ukb.pheno_to_mdd %>% 
  .[.$exposure %in% ls.sig.pheno_to_mdd$exposure,] %>% 
  as.data.frame %>% 
  select(-id.exposure,-id.outcome) %>% 
  mutate(Field=label_exposure) 


### Find traits using IVW only

mr.ivw.sig.ukb.pheno_to_mdd = mr.ukb.pheno_to_mdd %>% 
  filter(method=='Inverse variance weighted',p.fdr<0.05)


save(mr.sig.ukb.pheno_to_mdd,mr.ivw.sig.ukb.pheno_to_mdd,file='data/MR/traits_for_mvmr.RData')

write.table(mr.ivw.sig.ukb.pheno_to_mdd$exposure,file='data/MR/mvmr_ivwT_pheno.tsv',sep='\n',quote=F,row.names=F,col.names=F)
write.table(mr.sig.ukb.pheno_to_mdd$exposure,file='data/MR/mvmr_3methodsT_pheno.tsv',sep='\n',quote=F,row.names=F,col.names=F)



# Find SNPs for mvmr ------------------------------------------------------

ls.3method = mr.sig.ukb.pheno_to_mdd$exposure %>% paste0(.,'_TO_MDD')
ls.ivw = mr.ivw.sig.ukb.pheno_to_mdd$exposure %>% paste0(.,'_TO_MDD')

f.harmonisedDat.ivw = list.files(path='results/phewas/MR/ukb',pattern = '_mrDat.csv',full.names = T) %>% 
  .[grepl(paste0(ls.ivw,collapse = '|'),.)]

f.harmonisedDat.3method = list.files(path='results/phewas/MR/ukb',pattern = '_mrDat.csv',full.names = T) %>% 
  .[grepl(paste0(ls.3method,collapse = '|'),.)]

snps.ivw = f.harmonisedDat.ivw %>% as.list %>% 
  lapply(read_tsv) %>% 
  bind_rows %>% 
  {unique(.$SNP)}

snps.3methods = f.harmonisedDat.3method %>% as.list %>% 
  lapply(read_tsv) %>% 
  bind_rows %>% 
  {unique(.$SNP)}


write.table(snps.ivw,file='data/MR/mvmr_ivwT_snp.tsv',sep='\n',quote=F,row.names=F,col.names=F)
write.table(snps.3methods,file='data/MR/mvmr_3methodsT_snp.tsv',sep='\n',quote=F,row.names=F,col.names=F)
