PRS and MR PheWAS results
================
X Shen
06 January, 2022

-----

### PRS PheWAS

  - [Summary
    report](https://github.com/psychiatric-genomics-consortium/mdd-meta/blob/phewas/scripts/phewas/SUMMARY/summary.prs_phewas.md)

  - Result tables separated by category:
    ***results/phewas/phewas\_out\_\****

  - Models: ***results/phewas/models.rds***

  - Phenotypes used in the PheWAS:
    ***results/phewas/data\_dictionary/fields.final\****
    
    Use the following commands to retrieve a complete list of
    phenotypes:
    
    ``` r
    f.dic_dir = 'results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt'
    
    d.path = f.dic_dir %>% strsplit(.,'/') %>% unlist %>% 
      .[1:(length(.)-1)] %>% 
      paste0(.,collapse = '/')
    
    fields.all = list.files(d.path,pattern='^fields.final',full.names = T,include.dirs = T) %>% 
      as.list %>%
      lapply(fread,stringsAsFactors=F,header=T,sep='\t',quote="") %>%
      lapply(as.data.frame) %>% 
      lapply(mutate,field_used=as.character(field_used),FieldID=as.character(FieldID)) %>% 
      Reduce(function(dtf1,dtf2) bind_rows(dtf1,dtf2[,colnames(dtf1)]), .)
    ```

-----

### MR PheWAS

#### Analyses

  - Main analyses were conducted using clumped genome-wide significant
    genetic instruments for MDD and other phenotypes (p \< 5e-8).
    Bidiretional MR was conducted per MDD-phenotype pair. Analysis was
    first conducted using UKB phenotypes. Supplementary analysis used
    GWAS sumstats from [*the MRC IEU OpenGWAS
    database*](https://gwas.mrcieu.ac.uk/).

  - Sensitivity analysis 1: only the top 121 genetic instruments for MDD
    were used (instead of the ~500 genome-wide significant hits.) The
    number of genetic instruments was selected because it was the median
    number of the genetic instruments available for other phenotypes.

  - Sensitivity analysis 2: clumped genome-wide significant genetic
    instruments were used from the Howard et al MDD GWAS sumstats and
    other phenotypes (p \< 5e-8).

#### Summary reports and result files

  - Summary reports for the [main
    analyses](https://github.com/psychiatric-genomics-consortium/mdd-meta/blob/phewas/scripts/phewas/SUMMARY/summary.mr_phewas.md),
    [sensitivity
    analysis](https://github.com/psychiatric-genomics-consortium/mdd-meta/blob/phewas/scripts/phewas/SUMMARY/Sensitivity_analy.mr_phewas.md)
    and [scatter plots for
    QC](https://github.com/psychiatric-genomics-consortium/mdd-meta/blob/phewas/scripts/phewas/SUMMARY/QC.scatterplot.mr_phewas.md).

  - MR result tables: ***results/phewas/MR/\*/ALL\_mr\_res.tsv***

  - Data used for each individual analysis:
    ***results/phewas/MR/\*/\*\_mrDat.csv***

  - Individual scatter plot (adjustable) for each test:
    ***results/phewas/MR/\*/\*\_scatterplot.RData***

  - Individual scatter, funnel and leave-one-out plots for each test:
    ***results/phewas/MR/\*/\*\_plot.png***