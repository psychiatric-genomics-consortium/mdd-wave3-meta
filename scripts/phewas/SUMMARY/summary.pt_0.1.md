Phewas result
================
X Shen
11 August, 2021

    ## 
    ## Attaching package: 'kableExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

-----

## Load data

Inputs:

``` r
f.dat_dir = '/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Collab/mdd-meta/results/phewas/phewas_out_Body_MRI.rds'
f.dic_dir = '/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Collab/mdd-meta/results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt'
target.pT = '0.1'
f.category = '/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Collab/mdd-meta/data/phewas_categories.tsv'
```

Load
    results:

    ## Warning in FUN(X[[i]], ...): Found and resolved improper quoting out-of-
    ## sample. First healed line 680: <<"Assessment Centre > Imaging > Brain MRI >
    ## T1 structural brain MRI" 110 25756 "Scanner lateral (X) brain position" 43801
    ## 47219 "Accruing" "Continuous" "Units" "Data" "Auxiliary" "Unisex" 2 1 1317 "The
    ## exact location of the head and the radio-frequency receive coil in the scanner
    ## can affect data quality and imaging-derived phenotypes. To help account for
    ## variations in position in different scanned participants, several variables have
    ## been generated which describe aspects of the positioning and ca>>. If the fields
    ## are not quoted (e.g. field separator does not appear within any field), try
    ## quote="" to avoid this warning.

    ## Warning in FUN(X[[i]], ...): Found and resolved improper quoting in first 100
    ## rows. If the fields are not quoted (e.g. field separator does not appear within
    ## any field), try quote="" to avoid this warning.

### Summary

At pT=0.1 - Significant associations with the MDD3 PRS: 320 -
Significant associations with the Howard PRS143
