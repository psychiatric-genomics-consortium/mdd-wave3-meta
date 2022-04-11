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
library(tidyr)
library(ggplot2)
```

Read in LDSC rg tables

``` r
ldsc_rg_info <- read_tsv(snakemake@input$full) %>%
    mutate(dataset=str_match(p1, 'pgc_mdd_([:alpha:]+)_')[,2])
```

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   p1 = col_character(),
    ##   p2 = col_character(),
    ##   id = col_character(),
    ##   trait = col_character(),
    ##   subcategory = col_character(),
    ##   note = col_character(),
    ##   category = col_character(),
    ##   consortium = col_character(),
    ##   author = col_character(),
    ##   unit = col_character(),
    ##   population = col_character(),
    ##   sex = col_character(),
    ##   ontology = col_character(),
    ##   build = col_character(),
    ##   group_name = col_character(),
    ##   doi = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
ldsc_rg_mr_candidates <- read_tsv(snakemake@input$mr)
```

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   id = col_character(),
    ##   trait = col_character(),
    ##   rg = col_double(),
    ##   p = col_double(),
    ##   qvalue = col_double(),
    ##   gcov_int = col_double(),
    ##   subcategory = col_character()
    ## )

Plot size of genetic correlation versus covariance intercept

``` r
ggplot(ldsc_rg_info %>% filter(dataset=='full'), aes(x=abs(rg), y=gcov_int)) +
geom_point()
```

![](/Users/mark/Work/mdd-meta/docs/ldsc_files/figure-gfm/ldsc_rg_gcov-1.png)<!-- -->

Compare genetic covariance for `full` versus `noUKBB` sumstats

``` r
ldsc_full_vs_noukbb_gcov <-
ldsc_rg_info %>%
select(id, trait, dataset, gcov_int) %>%
pivot_wider(names_from=dataset, values_from=gcov_int) %>%
mutate(ukb=str_detect(id, 'ukb'))

ggplot(ldsc_full_vs_noukbb_gcov, aes(x=full, y=noUKBB, colour=ukb)) +
geom_point() + 
facet_grid(~ukb) +
coord_equal()
```

    ## Warning: Removed 7 rows containing missing values (geom_point).

![](/Users/mark/Work/mdd-meta/docs/ldsc_files/figure-gfm/ldsc_noukbb_gcov_point-1.png)<!-- -->

``` r
ldsc_full_noukbb_gcov <-
ldsc_rg_info %>%
filter(dataset %in% c('full', 'noUKBB')) %>%
mutate(ukb=if_else(str_detect(id, 'ukb'), 'sumstats from UKBB', 'other sumstats'))

ggplot(ldsc_full_noukbb_gcov, aes(x=gcov_int, fill=dataset)) +
geom_density() +
facet_grid(dataset~ukb)
```

![](/Users/mark/Work/mdd-meta/docs/ldsc_files/figure-gfm/ldsc_noukbb_gcov_hist-1.png)<!-- -->

``` r
ldsc_full_noukbb_gcov %>%
group_by(dataset, ukb) %>%
summarize(max_gov=max(abs(gcov_int)))
```

    ## `summarise()` has grouped output by 'dataset'. You can override using the `.groups` argument.

    ## # A tibble: 4 × 3
    ## # Groups:   dataset [2]
    ##   dataset ukb                max_gov
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

    ## # A tibble: 20 × 6
    ## # Groups:   subcategory [20]
    ##    trait                                  rg         p    qvalue gcov_int subcategory
    ##    <chr>                               <dbl>     <dbl>     <dbl>    <dbl> <chr>      
    ##  1 Neuroticism                        0.702  2.02e-162 4.89e-161   0.0446 <NA>       
    ##  2 Neuroticism                        0.658  3.36e-194 1.63e-192   0.0446 Personality
    ##  3 bipolar disorder                   0.402  1.10e- 77 9.65e- 77   0.0468 Psychiatri…
    ##  4 smoking initiation                 0.341  9.66e- 92 1.04e- 90   0.0479 Behavioural
    ##  5 Waist-to-hip ratio                 0.181  5.11e- 10 8.43e- 10  -0.0056 Anthropome…
    ##  6 Transferrin                        0.176  1.77e-  2 1.03e-  2  -0.017  Metal      
    ##  7 triglycerides                      0.159  2.85e- 17 8.12e- 17   0.0111 Lipid      
    ##  8 Years of schooling                -0.158  1.96e- 24 9.49e- 24  -0.0219 Education  
    ##  9 telomere length                   -0.124  5.33e- 10 8.75e- 10  -0.0033 Aging      
    ## 10 Age at menarche                   -0.117  3.42e- 10 5.81e- 10  -0.0082 Reproducti…
    ## 11 C-Reactive protein level           0.106  5.26e-  6 6.33e-  6   0.0083 Immune sys…
    ## 12 ER+ Breast cancer (GWAS)           0.103  1   e-  4 9.59e-  5  -0.01   Cancer     
    ## 13 Crohn's disease                    0.102  3.27e-  6 4.05e-  6  -0.0027 Autoimmune…
    ## 14 Urinary sodium-potassium ratio     0.0817 1.1 e-  3 8.95e-  4  -0.0104 Biomarker  
    ## 15 Heart rate                         0.0767 1.1 e-  2 6.91e-  3  -0.0078 Hemodynamic
    ## 16 Platelet count                     0.068  9.6 e-  3 6.14e-  3  -0.0004 Haemotolog…
    ## 17 Fasting glucose                    0.0607 4.85e-  2 2.63e-  2  -0.0068 Glycemic   
    ## 18 Chronotype                        -0.0594 1.34e-  2 8.17e-  3  -0.0107 Sleeping   
    ## 19 Femoral neck bone mineral density -0.0578 3.07e-  2 1.72e-  2   0.0026 Bone       
    ## 20 diastolic blood pressure           0.0348 9.3 e-  3 5.97e-  3  -0.0185 Blood pres…

Plot rg for known subcategories

``` r
ggplot(ldsc_rg_info %>% filter(!is.na(subcategory), dataset %in% 'full'),
aes(x=trait, y=rg, ymin=rg+se*qnorm(0.025), ymax=rg+se*qnorm(0.975))) +
geom_pointrange() +
coord_flip()
```

![](/Users/mark/Work/mdd-meta/docs/ldsc_files/figure-gfm/rg_subcats-1.png)<!-- -->

Compare new and previous results

``` r
ldsc_rg_full_howard <-
ldsc_rg_info %>%
filter(dataset %in% c('full', 'howard')) %>%
select(dataset, rg, se, p, gcov_int, id, trait) %>%
pivot_wider(id_cols=c(id, trait),
            names_from=dataset,
            values_from=c(rg, se, p, gcov_int))
```

Results that were `NA` with the previous sumstats

``` r
ldsc_rg_full_howard %>%
filter(!is.na(p_full), is.na(p_howard))
```

    ## # A tibble: 1 × 10
    ##   id     trait rg_full rg_howard se_full se_howard p_full p_howard gcov_int_full
    ##   <chr>  <chr>   <dbl>     <dbl>   <dbl>     <dbl>  <dbl>    <dbl>         <dbl>
    ## 1 ebi-a… CD27…   0.106        NA   0.170        NA  0.534       NA        0.0048
    ## # … with 1 more variable: gcov_int_howard <dbl>

Traits that were not significantly correlated before, after FDR
correction

``` r
ldsc_rg_full_howard_fdr <-
ldsc_rg_full_howard %>%
filter(!is.na(p_full), !is.na(p_howard)) %>%
mutate(q_full=fdrtool::fdrtool(p_full, statistic='p', plot=FALSE)$qval,
       q_howard=fdrtool::fdrtool(p_howard, statistic='p', plot=FALSE)$qval)
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr
    ## 
    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

``` r
ldsc_rg_full_howard_fdr %>%
filter(p_full <= 0.05, q_full <= 0.05, q_howard > 0.05) %>%
select(trait, rg_full, q_full) %>%
arrange(desc(abs(rg_full))) %>%
print(n=Inf)
```

    ## # A tibble: 21 × 3
    ##    trait                                                        rg_full   q_full
    ##    <chr>                                                          <dbl>    <dbl>
    ##  1 Added milk to instant coffee                                 -0.165   6.52e-3
    ##  2 Added milk to standard tea                                   -0.112   3.26e-3
    ##  3 Alcohol drinker status: Never                                -0.0888  3.65e-3
    ##  4 Pulse rate (during blood-pressure measurement)                0.0857  2.20e-3
    ##  5 Relative age voice broke                                      0.0819  4.52e-4
    ##  6 HDL cholesterol                                              -0.0771  1.81e-3
    ##  7 Corneal resistance factor (left)                              0.0705  7.24e-5
    ##  8 Mineral and other dietary supplements: Glucosamine           -0.0589  5.14e-3
    ##  9 Mineral and other dietary supplements: Fish oil (including …  0.0563  5.02e-3
    ## 10 Types of transport used (excluding work): Walk               -0.0527  2.46e-3
    ## 11 Relative age of first facial hair                            -0.0486  4.54e-3
    ## 12 Waist circumference                                           0.0445  4.77e-3
    ## 13 Waist circumference                                           0.0432  5.46e-3
    ## 14 Mean corpuscular hemoglobin                                   0.0425  1.13e-3
    ## 15 Mean reticulocyte volume                                      0.0404  6.41e-4
    ## 16 Red cell distribution width                                   0.0397  8.18e-3
    ## 17 Red blood cell (erythrocyte) distribution width               0.0382  2.62e-3
    ## 18 Mean corpuscular volume                                       0.0377  2.85e-3
    ## 19 Heel bone mineral density (BMD) T-score  automated (left)     0.0355  8.82e-3
    ## 20 diastolic blood pressure                                      0.0348  1.98e-3
    ## 21 Mean corpuscular volume                                       0.0304  4.70e-3

Non-UKB/overlapping phenotypes sorted by magnitude of genetic
correlation

``` r
ldsc_rg_full_howard_fdr %>%
    filter(!str_detect(id, 'ukb'), p_full <= 0.05, q_full <= 0.05, abs(gcov_int_full) <= 0.05) %>%
    arrange(desc(abs(rg_full))) %>%
    select(trait, rg_full, id, p=p_full, q=q_full) %>%
    print(n=40)
```

    ## # A tibble: 177 × 5
    ##    trait                                  rg_full id                 p         q
    ##    <chr>                                    <dbl> <chr>          <dbl>     <dbl>
    ##  1 Neuroticism                              0.702 ebi-a-GCS… 2.02e-162 5.44e-162
    ##  2 Neuroticism                              0.671 ebi-a-GCS… 8.31e-164 2.28e-163
    ##  3 Neuroticism                              0.658 ieu-a-1007 3.36e-194 1.20e-193
    ##  4 Wellbeing                               -0.631 ieu-b-4852 1.27e- 26 7.04e- 27
    ##  5 Feeling tense                            0.585 ebi-a-GCS… 8.55e-219 3.80e-218
    ##  6 Irritable mood                           0.496 ebi-a-GCS… 6.67e-143 1.50e-142
    ##  7 Feeling worry                            0.447 ebi-a-GCS… 1.96e-144 4.48e-144
    ##  8 Feeling nervous                          0.403 ebi-a-GCS… 2.25e- 85 3.02e- 85
    ##  9 bipolar disorder                         0.402 ieu-b-41   1.10e- 77 1.37e- 77
    ## 10 Trauma exposure in major depressive d…   0.400 ebi-a-GCS… 3.16e- 24 1.66e- 24
    ## 11 Ever smoked                              0.390 ieu-b-4858 1.16e- 50 1.07e- 50
    ## 12 schizophrenia                            0.374 ieu-b-42   1.57e- 97 2.36e- 97
    ## 13 Schizophrenia                            0.365 ieu-a-22   7.32e- 98 1.13e- 97
    ## 14 Bipolar disorder                         0.356 ieu-a-801  3.96e- 24 2.08e- 24
    ## 15 Worry too long after an embarrassing …   0.343 ebi-a-GCS… 3.17e- 63 3.39e- 63
    ## 16 smoking initiation                       0.341 ieu-b-4877 9.66e- 92 1.39e- 91
    ## 17 Age Of Smoking Initiation               -0.307 ieu-b-24   6.73e- 41 5.12e- 41
    ## 18 Cigarettes per Day                       0.303 ieu-b-25   1.83e- 23 9.46e- 24
    ## 19 Cigarettes smoked per day                0.303 ieu-b-142  1.83e- 23 9.44e- 24
    ## 20 Parental longevity (father's age at d…  -0.285 ebi-a-GCS… 2.42e- 21 1.17e- 21
    ## 21 Coronary artery disease                  0.278 ebi-a-GCS… 6.98e- 40 5.20e- 40
    ## 22 Strenuous sports or other exercises     -0.274 ebi-a-GCS… 3.79e- 40 2.84e- 40
    ## 23 Coronary artery disease                  0.246 ebi-a-GCS… 2.58e- 45 2.13e- 45
    ## 24 Type 2 diabetes                          0.215 ebi-a-GCS… 1.38e- 24 7.40e- 25
    ## 25 Anorexia Nervosa                         0.213 ieu-a-1186 3.03e-  7 9.13e-  8
    ## 26 C-reactive protein                       0.211 ieu-b-4764 1   e-  4 2.49e-  5
    ## 27 Knee osteoarthritis                      0.206 ebi-a-GCS… 1.61e- 20 7.67e- 21
    ## 28 Parental longevity (combined parental…   0.205 ebi-a-GCS… 4.52e- 14 1.82e- 14
    ## 29 25 hydroxyvitamin D level               -0.204 ieu-b-4808 1.40e- 22 7.04e- 23
    ## 30 College completion                      -0.200 ieu-a-836  2.09e- 13 8.23e- 14
    ## 31 Allergic disease (asthma, hay fever o…   0.197 ebi-a-GCS… 6.03e- 19 2.77e- 19
    ## 32 Years of schooling                      -0.196 ieu-a-755  4.43e- 13 1.72e- 13
    ## 33 Years of schooling                      -0.184 ieu-a-1010 8.99e- 18 3.96e- 18
    ## 34 25 hydroxyvitamin D level               -0.181 ieu-b-4812 2.91e- 19 1.34e- 19
    ## 35 Waist-to-hip ratio                       0.181 ieu-a-77   5.11e- 10 1.74e- 10
    ## 36 Waist-to-hip ratio                       0.181 ieu-a-76   2.29e- 10 7.94e- 11
    ## 37 Years of schooling                      -0.179 ieu-a-1001 1.08e- 22 5.45e- 23
    ## 38 Transferrin                              0.176 ieu-a-1052 1.77e-  2 3.63e-  3
    ## 39 Years of schooling                      -0.173 ieu-b-4836 3.10e- 10 1.07e- 10
    ## 40 Accelerometer-based physical activity…  -0.172 ebi-a-GCS… 5.61e- 13 2.16e- 13
    ## # … with 137 more rows
