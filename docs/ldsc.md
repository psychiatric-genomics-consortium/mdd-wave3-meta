LD Score Genetic Correlations
================

``` r
library(readr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(stringr)
library(ggplot2)
```

Read in LDSC rg tables

``` r
ldsc_rg_info <- read_tsv(snakemake@input$full)
```

    ## Rows: 1074 Columns: 35

    ## ── Column specification ─────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): p1, p2, id, trait, subcategory, note, category, consortium, author...
    ## dbl (19): rg, se, z, p, h2_obs, h2_obs_se, h2_int, h2_int_se, gcov_int, gcov...

    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ldsc_noukbb_rg_info <- read_tsv(snakemake@input$noukbb)
```

    ## Rows: 1075 Columns: 35

    ## ── Column specification ─────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): p1, p2, id, trait, subcategory, note, category, consortium, author...
    ## dbl (19): rg, se, z, p, h2_obs, h2_obs_se, h2_int, h2_int_se, gcov_int, gcov...

    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ldsc_rg_mr_candidates <- read_tsv(snakemake@input$mr)
```

    ## Rows: 172 Columns: 7

    ## ── Column specification ─────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): id, trait, subcategory
    ## dbl (4): rg, p, qvalue, gcov_int

    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Plot size of genetic correlation versus covariance intercept

``` r
ggplot(ldsc_rg_info, aes(x=abs(rg), y=gcov_int)) +
geom_point()
```

![](/gpfs/igmmfs01/eddie/GenScotDepression/madams23/projects/mdd-meta/docs/ldsc_files/figure-gfm/ldsc_rg_gcov-1.png)<!-- -->

Compare genetic covariance for `full` versus `noUKBB` sumstats

``` r
ldsc_full_vs_noukbb_gcov <-
inner_join(
ldsc_rg_info %>% transmute(id, trait, gcov_int.full=gcov_int),
ldsc_noukbb_rg_info %>% transmute(id, trait, gcov_int.noUKBB=gcov_int),
by=c('id', 'trait')
) %>%
mutate(ukb=str_detect(id, 'ukb'))

ggplot(ldsc_full_vs_noukbb_gcov, aes(x=gcov_int.full, y=gcov_int.noUKBB, colour=ukb)) +
geom_point() + 
facet_grid(~ukb) +
coord_equal()
```

![](/gpfs/igmmfs01/eddie/GenScotDepression/madams23/projects/mdd-meta/docs/ldsc_files/figure-gfm/ldsc_noukbb_gcov_point-1.png)<!-- -->

``` r
ldsc_full_noukbb_gcov <-
bind_rows(
ldsc_rg_info %>% transmute(id, trait, gcov_int, cohorts='full'),
ldsc_noukbb_rg_info %>% transmute(id, trait, gcov_int=gcov_int, cohorts='noUKBB')
) %>%
mutate(ukb=if_else(str_detect(id, 'ukb'), 'sumstats from UKBB', 'other sumstats'))

ggplot(ldsc_full_noukbb_gcov, aes(x=gcov_int, fill=cohorts)) +
geom_density() +
facet_grid(cohorts~ukb)
```

![](/gpfs/igmmfs01/eddie/GenScotDepression/madams23/projects/mdd-meta/docs/ldsc_files/figure-gfm/ldsc_noukbb_gcov_hist-1.png)<!-- -->

``` r
ldsc_full_noukbb_gcov %>%
group_by(cohorts, ukb) %>%
summarize(max_gov=max(abs(gcov_int)))
```

    ## `summarise()` has grouped output by 'cohorts'. You can override using the `.groups` argument.

    ## # A tibble: 4 × 3
    ## # Groups:   cohorts [2]
    ##   cohorts ukb                max_gov
    ##   <chr>   <chr>                <dbl>
    ## 1 full    other sumstats       0.266
    ## 2 full    sumstats from UKBB   0.176
    ## 3 noUKBB  other sumstats       0.245
    ## 4 noUKBB  sumstats from UKBB   0.032

FDR corrected associations, removing phenotypes from UKB and with large
genetic covariance intercepts

``` r
ldsc_rg_mr_candidates %>%
group_by(subcategory) %>%
filter(rg==max(rg)) %>%
filter(p==min(p)) %>%
select(-id)
```

    ## # A tibble: 19 × 6
    ## # Groups:   subcategory [19]
    ##    trait                               rg        p   qvalue gcov_int subcategory
    ##    <chr>                            <dbl>    <dbl>    <dbl>    <dbl> <chr>      
    ##  1 Trauma exposure in major depr…  0.400  3.16e-24 4.37e-23  -0.0011 <NA>       
    ##  2 Cigarettes smoked per day       0.303  1.83e-23 1.97e-22   0.0237 Behavioural
    ##  3 Anorexia Nervosa                0.213  3.03e- 7 5.56e- 7   0.0118 Psychiatri…
    ##  4 Waist-to-hip ratio              0.181  5.11e-10 1.16e- 9  -0.0056 Anthropome…
    ##  5 Transferrin                     0.176  1.77e- 2 1.15e- 2  -0.017  Metal      
    ##  6 triglycerides                   0.159  2.85e-17 1.45e-16   0.0111 Lipid      
    ##  7 Years of schooling             -0.158  1.96e-24 3.04e-23  -0.0219 Education  
    ##  8 telomere length                -0.124  5.33e-10 1.20e- 9  -0.0033 Aging      
    ##  9 Age at menarche                -0.117  3.42e-10 8.08e-10  -0.0082 Reproducti…
    ## 10 C-Reactive protein level        0.106  5.26e- 6 7.90e- 6   0.0083 Immune sys…
    ## 11 ER+ Breast cancer (GWAS)        0.103  1   e- 4 1.14e- 4  -0.01   Cancer     
    ## 12 Crohn's disease                 0.102  3.27e- 6 5.10e- 6  -0.0027 Autoimmune…
    ## 13 Urinary sodium-potassium ratio  0.0817 1.1 e- 3 1.03e- 3  -0.0104 Biomarker  
    ## 14 Heart rate                      0.0767 1.1 e- 2 7.77e- 3  -0.0078 Hemodynamic
    ## 15 Platelet count                  0.068  9.6 e- 3 6.92e- 3  -0.0004 Haemotolog…
    ## 16 Fasting glucose                 0.0607 4.85e- 2 2.91e- 2  -0.0068 Glycemic   
    ## 17 Chronotype                     -0.0594 1.34e- 2 9.15e- 3  -0.0107 Sleeping   
    ## 18 Femoral neck bone mineral den… -0.0578 3.07e- 2 1.91e- 2   0.0026 Bone       
    ## 19 diastolic blood pressure        0.0348 9.3 e- 3 6.73e- 3  -0.0185 Blood pres…

Plot rg for known subcategories

``` r
ggplot(ldsc_rg_info %>% filter(!is.na(subcategory)), aes(x=trait, y=rg, ymin=rg+se*qnorm(0.025), ymax=rg+se*qnorm(0.975))) +
geom_pointrange() +
coord_flip()
```

![](/gpfs/igmmfs01/eddie/GenScotDepression/madams23/projects/mdd-meta/docs/ldsc_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
