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
    ## ── Column specification ─────────────────────────────────────────────────────────────────────
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
    ## ── Column specification ─────────────────────────────────────────────────────────────────────
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
select(dataset, rg, se, p, gcov_int, id, trait, subcategory) %>%
pivot_wider(id_cols=c(subcategory, id, trait),
            names_from=dataset,
            values_from=c(rg, se, p, gcov_int))
```

Results that were `NA` with the previous sumstats

``` r
ldsc_rg_full_howard %>%
filter(!is.na(p_full), is.na(p_howard))
```

    ## # A tibble: 1 × 11
    ##   subcategory id      trait  rg_full rg_howard se_full se_howard p_full p_howard
    ##   <chr>       <chr>   <chr>    <dbl>     <dbl>   <dbl>     <dbl>  <dbl>    <dbl>
    ## 1 <NA>        ebi-a-… CD27 …   0.106        NA   0.170        NA  0.534       NA
    ## # … with 2 more variables: gcov_int_full <dbl>, gcov_int_howard <dbl>

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
select(trait, rg_full, p_full, q_full) %>%
arrange(desc(abs(rg_full))) %>%
print(n=Inf)
```

    ## # A tibble: 21 × 4
    ##    trait                                                 rg_full p_full   q_full
    ##    <chr>                                                   <dbl>  <dbl>    <dbl>
    ##  1 Added milk to instant coffee                          -0.165  0.0327  6.52e-3
    ##  2 Added milk to standard tea                            -0.112  0.0158  3.26e-3
    ##  3 Alcohol drinker status: Never                         -0.0888 0.0178  3.65e-3
    ##  4 Pulse rate (during blood-pressure measurement)         0.0857 0.0104  2.20e-3
    ##  5 Relative age voice broke                               0.0819 0.002   4.52e-4
    ##  6 HDL cholesterol                                       -0.0771 0.0085  1.81e-3
    ##  7 Corneal resistance factor (left)                       0.0705 0.0003  7.24e-5
    ##  8 Mineral and other dietary supplements: Glucosamine    -0.0589 0.0255  5.14e-3
    ##  9 Mineral and other dietary supplements: Fish oil (inc…  0.0563 0.0249  5.02e-3
    ## 10 Types of transport used (excluding work): Walk        -0.0527 0.0117  2.46e-3
    ## 11 Relative age of first facial hair                     -0.0486 0.0224  4.54e-3
    ## 12 Waist circumference                                    0.0445 0.0236  4.77e-3
    ## 13 Waist circumference                                    0.0432 0.0272  5.46e-3
    ## 14 Mean corpuscular hemoglobin                            0.0425 0.0052  1.13e-3
    ## 15 Mean reticulocyte volume                               0.0404 0.0029  6.41e-4
    ## 16 Red cell distribution width                            0.0397 0.0414  8.18e-3
    ## 17 Red blood cell (erythrocyte) distribution width        0.0382 0.0125  2.62e-3
    ## 18 Mean corpuscular volume                                0.0377 0.0137  2.85e-3
    ## 19 Heel bone mineral density (BMD) T-score  automated (…  0.0355 0.0448  8.82e-3
    ## 20 diastolic blood pressure                               0.0348 0.0093  1.98e-3
    ## 21 Mean corpuscular volume                                0.0304 0.0232  4.70e-3

Non-UKB/overlapping phenotypes sorted by magnitude of genetic
correlation

``` r
ldsc_rg_full_howard_fdr %>%
    filter(!str_detect(id, 'ukb'), p_full <= 0.05, q_full <= 0.05, abs(gcov_int_full) <= 0.05) %>%
    group_by(subcategory, word(str_to_lower(trait), 1)) %>%
    filter(p_full == min(p_full)) %>% 
    filter(rg_full == max(rg_full)) %>%
    arrange(subcategory, desc(abs(rg_full))) %>%
    group_by(subcategory) %>%
    slice(1:6) %>%
    select(trait, rg=rg_full, se=se_full, FDR=q_full, id) %>%
    print(n=Inf)
```

    ## Adding missing grouping variables: `subcategory`

    ## # A tibble: 42 × 6
    ## # Groups:   subcategory [20]
    ##    subcategory                trait                  rg     se       FDR id     
    ##    <chr>                      <chr>               <dbl>  <dbl>     <dbl> <chr>  
    ##  1 Aging                      telomere length   -0.124  0.0199 1.81e- 10 ieu-b-…
    ##  2 Anthropometric             Waist-to-hip rat…  0.181  0.0285 7.94e- 11 ieu-a-…
    ##  3 Anthropometric             Hip circumference  0.142  0.0203 9.76e- 13 ieu-a-…
    ##  4 Anthropometric             Waist circumfere…  0.142  0.0202 7.91e- 13 ieu-a-…
    ##  5 Anthropometric             body mass index    0.141  0.0141 1.41e- 23 ieu-b-…
    ##  6 Anthropometric             Obesity class 1    0.105  0.0204 7.33e-  8 ieu-a-…
    ##  7 Anthropometric             Height            -0.0518 0.0158 2.33e-  4 ieu-a-…
    ##  8 Autoimmune / inflammatory  Crohn's disease    0.102  0.0218 9.04e-  7 ieu-a-…
    ##  9 Behavioural                smoking initiati…  0.341  0.0168 1.39e- 91 ieu-b-…
    ## 10 Behavioural                Age Of Smoking I… -0.307  0.0229 5.12e- 41 ieu-b-…
    ## 11 Behavioural                Cigarettes smoke…  0.303  0.0304 9.44e- 24 ieu-b-…
    ## 12 Behavioural                Alcoholic drinks…  0.079  0.0199 1.77e-  5 ieu-b-…
    ## 13 Biomarker                  Urinary sodium-p…  0.0817 0.025  2.55e-  4 ieu-b-…
    ## 14 Blood pressure             diastolic blood …  0.0348 0.0134 1.98e-  3 ieu-b-…
    ## 15 Bone                       Femoral neck bon… -0.0578 0.0267 6.14e-  3 ieu-a-…
    ## 16 Cancer                     ER+ Breast cance…  0.103  0.0266 2.49e-  5 ieu-a-…
    ## 17 Cancer                     Breast cancer (i…  0.0936 0.0254 4.88e-  5 ieu-a-…
    ## 18 Education                  College completi… -0.200  0.0273 8.23e- 14 ieu-a-…
    ## 19 Education                  Years of schooli… -0.158  0.0155 1.04e- 24 ieu-a-…
    ## 20 Glycemic                   Fasting glucose    0.0607 0.0308 9.52e-  3 ieu-b-…
    ## 21 Haemotological             Platelet count     0.068  0.0262 2.04e-  3 ieu-a-…
    ## 22 Hemodynamic                Heart rate         0.0767 0.0302 2.32e-  3 ieu-a-…
    ## 23 Immune system              C-Reactive prote…  0.106  0.0232 1.43e-  6 ieu-b-…
    ## 24 Lipid                      triglycerides      0.159  0.0188 1.23e- 17 ieu-b-…
    ## 25 Lipid                      HDL cholesterol   -0.130  0.0179 1.72e- 13 ieu-b-…
    ## 26 Lipid                      apolipoprotein A… -0.0953 0.0187 1.07e-  7 ieu-b-…
    ## 27 Lipid                      LDL cholesterol   -0.0674 0.0188 7.24e-  5 ieu-b-…
    ## 28 Metal                      Transferrin        0.176  0.0742 3.63e-  3 ieu-a-…
    ## 29 Personality                Neuroticism        0.658  0.0221 1.20e-193 ieu-a-…
    ## 30 Psychiatric / neurological bipolar disorder   0.402  0.0215 1.37e- 77 ieu-b-…
    ## 31 Psychiatric / neurological Schizophrenia      0.365  0.0174 1.13e- 97 ieu-a-…
    ## 32 Psychiatric / neurological Anorexia Nervosa   0.213  0.0416 9.13e-  8 ieu-a-…
    ## 33 Psychiatric / neurological multiple scleros…  0.110  0.0229 4.07e-  7 ieu-b-…
    ## 34 Reproductive aging         Age at menarche   -0.117  0.0186 1.18e- 10 ieu-a-…
    ## 35 Sleeping                   Sleep duration    -0.084  0.0297 1.00e-  3 ieu-a-…
    ## 36 Sleeping                   Chronotype        -0.0594 0.024  2.79e-  3 ieu-a-…
    ## 37 <NA>                       Neuroticism        0.671  0.0246 2.28e-163 ebi-a-…
    ## 38 <NA>                       Wellbeing         -0.631  0.0591 7.04e- 27 ieu-b-…
    ## 39 <NA>                       Feeling tense      0.585  0.0185 3.80e-218 ebi-a-…
    ## 40 <NA>                       Irritable mood     0.496  0.0195 1.50e-142 ebi-a-…
    ## 41 <NA>                       Trauma exposure …  0.400  0.0393 1.66e- 24 ebi-a-…
    ## 42 <NA>                       Ever smoked        0.390  0.0261 1.07e- 50 ieu-b-…
