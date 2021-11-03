PGC MDD3 sumstats QC checks for meta-analysis
================

-   [Sample sizes](#sample-sizes)
-   [Reference panel alignment](#reference-panel-alignment)
-   [Genetic correlation with clinical
    cohorts](#genetic-correlation-with-clinical-cohorts)
-   [Genetic correlation with previous meta
    analysis](#genetic-correlation-with-previous-meta-analysis)
-   [Heritabilities](#heritabilities)
-   [Genetic covariance intercepts](#genetic-covariance-intercepts)
    -   [Clustering](#clustering)
    -   [Genetic correlations](#genetic-correlations)

``` r
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
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(fdrtool)
library(corrplot)
```

    ## corrplot 0.90 loaded

# Sample sizes

``` r
meta_qc_align <- read_tsv(snakemake@input$meta_qc_align)
```

    ## 
    ## -- Column specification ----------------------------------------------------------------------------------------------------------
    ## cols(
    ##   .default = col_double(),
    ##   cohort = col_character(),
    ##   ancestries = col_character(),
    ##   release = col_character()
    ## )
    ## i Use `spec()` for the full column specifications.

``` r
cohorts_mdd <- read_tsv(snakemake@input$cohorts_mdd)
```

    ## Warning: Missing column names filled in: 'X7' [7]

    ## 
    ## -- Column specification ----------------------------------------------------------------------------------------------------------
    ## cols(
    ##   Dataset = col_character(),
    ##   N_cases = col_double(),
    ##   N_controls = col_double(),
    ##   `LAMBDA-GC` = col_double(),
    ##   `N-SNPs` = col_double(),
    ##   N_eff_half = col_double(),
    ##   X7 = col_logical()
    ## )

``` r
cohorts_samples <- cohorts_mdd %>%
filter(Dataset != 'SUM') %>%
mutate(cohort=str_match(Dataset, "mdd_(.+)_eur")[,2])

meta_samples <- 
meta_qc_align %>%
filter(!cohort %in% c('MDD29', 'PGC')) %>%
bind_rows(cohorts_samples) %>%
group_by(cohort) %>%
summarize(Cases=sum(N_cases), Controls=sum(N_controls)) %>%
pivot_longer(Cases:Controls, names_to="MDD", values_to="N") %>%
arrange(desc(N))

cohort_order <-
meta_samples %>%
filter(MDD == 'Cases') %>%
arrange(desc(N)) %>%
pull(cohort)

ggplot(meta_samples, aes(x=0, y=0, color=MDD, size=N)) +
geom_point() +
facet_wrap(~factor(cohort, levels=cohort_order)) +
scale_size_area(max_size=33) +
theme_minimal() +
theme(axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank())
```

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/meta_sample_size-1.png)<!-- -->

# Reference panel alignment

Marker names were aligned to the meta-analysis reference panel and
effect sizes were checked. SNPs with MAF \< 0.001 in the imputation
reference panel and INFO \< 0.1 in the cohort were removed.

``` r
meta_qc_align %>% print(n=Inf)
```

    ## # A tibble: 32 x 22
    ##    cohort   ancestries release            N_cases N_controls    snps snps_merged
    ##    <chr>    <chr>      <chr>                <dbl>      <dbl>   <dbl>       <dbl>
    ##  1 23andMe  eas        v7_2                  2727      90220  8.24e6     6538698
    ##  2 BBJ      eas        hum0197v3_Depv1_2~     836     177794  1.34e7     7685385
    ##  3 CONVERGE eas        10640                 5303       5337  4.64e6     4614299
    ##  4 Taiwan   eas        20200327              1348       8392  4.15e6     4153855
    ##  5 23andMe  eur        v7_2_202012         112892    1773938  1.33e7    10741251
    ##  6 AGDS     eur        202012               12123      12684  7.62e6     7583592
    ##  7 Airwave  eur        0820                  2100      15713  7.49e6     7473166
    ##  8 ALSPAC   eur        27022020               472       3475  9.09e6     7672087
    ##  9 BASIC    eur        202011                1003       1854  7.77e6     7731069
    ## 10 BioVU    eur        NoCov_SAIGE_051821    7757      24723  6.26e6     6260491
    ## 11 DBDS     eur        FINAL202103          13347     145996  7.59e6     6890112
    ## 12 deCODE   eur        DEPALL_FINAL_WHEAD   20000      28000  8.84e6     7763310
    ## 13 ESTBB    eur        EstBB                35473      91301  2.69e7    15112315
    ## 14 EXCEED   eur        202010                 580       2071  8.08e6     8046825
    ## 15 FinnGen  eur        R5_18032020          23424     192220  1.64e7    12958898
    ## 16 GenScot  eur        SCID_0721a             930       5730  7.80e6     7790538
    ## 17 GERA     eur        0915a_mds5            7162      38307  1.09e7     9020123
    ## 18 HUNT     eur        gp_hospital_metac~   11658      42535  8.69e6     7802064
    ## 19 iPSYCH   eur        2012_HRC             19156      22708  8.81e6     8773714
    ## 20 iPSYCH   eur        2015i_HRC            10002      15434  8.85e6     8814173
    ## 21 lgic2    eur        202011                 906       4717  7.76e6     7723129
    ## 22 MDD49    eur        29w2_20w3_X28w2_1~   28147      48033  7.94e6     7906555
    ## 23 MoBa     eur        harvest12              603      10213  6.50e6     6501052
    ## 24 MoBa     eur        harvest24              367       6122  6.50e6     6501041
    ## 25 MoBa     eur        rotterdam1             553       8860  6.50e6     6501057
    ## 26 MVP      eur        rel4icdDEP_Geno_2~  151974     226640  1.60e7    13720146
    ## 27 PBK      eur        2020                  5607      16080  7.35e6     7325429
    ## 28 PREFECT  eur        run1                  1796       3290  8.96e6     8919673
    ## 29 SHARE    eur        godartsshare_8420~    1063       1921  1.36e7    12504715
    ## 30 STAGE    eur        MDDdx_saige            421       9134  7.43e6     7398077
    ## 31 tkda1    eur        run1                   672        846  8.86e6     8817113
    ## 32 UKBB     eur        MD_glm_202107        54669     306461  1.19e7     9347557
    ## # ... with 15 more variables: snps_unambiguous_flips <dbl>,
    ## #   snps_matching <dbl>, snps_turned <dbl>, snps_unresolved <dbl>,
    ## #   snps_aligned <dbl>, median_fst <dbl>, max_fst <dbl>, var_fst <dbl>,
    ## #   snps_kept <dbl>, median_OR <dbl>, max_OR <dbl>, median_OR01 <dbl>,
    ## #   median_SE <dbl>, max_SE <dbl>, median_SE01 <dbl>

Median odds-ratios with standard error

``` r
ggplot(meta_qc_align %>%
       unite('cohort',
              cohort, ancestries, release, 
              sep="."),
        aes(x=reorder(cohort, (4*N_cases*N_controls)/(N_cases + N_controls)),
            y=median_OR01,
            ymax=exp(log(median_OR01)+median_SE01),
            ymin=exp(log(median_OR01)-median_SE01))) +
geom_linerange() +
geom_point() +
coord_flip(ylim=c(0.8, 1.25))
```

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/meta_qc_align_mean_OR-1.png)<!-- -->

Mean standard error versus effective sample size:

``` r
ggplot(meta_qc_align %>%
       unite('cohort',
              cohort, ancestries, release, 
              sep=".") %>%
       mutate(Neff=(4*N_cases*N_controls)/(N_cases + N_controls)),
        aes(x=Neff, y=median_SE01, label=cohort)) +
geom_text(size=2) +
scale_x_log10(breaks=c(1, 100, 1000, 5000, 10000, 50000, 100000, 500000))
```

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/meta_qc_align_Neff_SE-1.png)<!-- -->

# Genetic correlation with clinical cohorts

LDSC genetic correlations were calculated with the PGC `MDD29` clinical
cohorts

``` r
meta_qc_ldsc <-
read_table2(snakemake@input$meta_qc_ldsc) %>%
mutate(se.mdd2=as.numeric(str_remove_all(se.mdd2, "[\\(\\)]")),
       se.mdd29=as.numeric(str_remove_all(se.mdd29, "[\\(\\)]")),
       gencov_se.mdd2=as.numeric(str_remove_all(gencov_se.mdd2, "[\\(\\)]")),
       gencov_se.mdd29=as.numeric(str_remove_all(gencov_se.mdd29, "[\\(\\)]"))
    )
```

    ## 
    ## -- Column specification ----------------------------------------------------------------------------------------------------------
    ## cols(
    ##   cohort = col_character(),
    ##   release = col_character(),
    ##   h2_obs = col_double(),
    ##   h2_obs_se = col_double(),
    ##   rg.mdd2 = col_double(),
    ##   se.mdd2 = col_double(),
    ##   gencov.mdd2 = col_double(),
    ##   gencov_se.mdd2 = col_double(),
    ##   z1z2.mdd2 = col_double(),
    ##   rg.mdd29 = col_double(),
    ##   se.mdd29 = col_double(),
    ##   gencov.mdd29 = col_double(),
    ##   gencov_se.mdd29 = col_double(),
    ##   z1z2.mdd29 = col_double()
    ## )

    ## Warning: 3 parsing failures.
    ## row col   expected    actual                           file
    ##   4  -- 14 columns 7 columns 'docs/tables/meta_qc_ldsc.txt'
    ##   5  -- 14 columns 7 columns 'docs/tables/meta_qc_ldsc.txt'
    ##  20  -- 14 columns 7 columns 'docs/tables/meta_qc_ldsc.txt'

``` r
ggplot(meta_qc_ldsc, aes(x=reorder(paste(cohort, release), rg.mdd29), y=rg.mdd29, ymin=rg.mdd29-se.mdd29, ymax=rg.mdd29+se.mdd29)) +
geom_linerange() +
geom_point() +
coord_flip(ylim=c(0, 1.25))
```

    ## Warning: Removed 8 rows containing missing values (geom_segment).

    ## Warning: Removed 8 rows containing missing values (geom_point).

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/metq_qc_ldsc_rg-1.png)<!-- -->

Genetic covariance:

``` r
ggplot(meta_qc_ldsc, aes(x=reorder(paste(cohort, release), gencov.mdd29), y=gencov.mdd29, ymin=gencov.mdd29-gencov_se.mdd29, ymax=gencov.mdd29+gencov_se.mdd29)) +
geom_linerange() +
geom_point() +
coord_flip()
```

    ## Warning: Removed 3 rows containing missing values (geom_segment).

    ## Warning: Removed 3 rows containing missing values (geom_point).

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/metq_qc_ldsc_gcov-1.png)<!-- -->

# Genetic correlation with previous meta analysis

``` r
ggplot(meta_qc_ldsc, aes(x=reorder(paste(cohort, release), rg.mdd2), y=rg.mdd2, ymin=rg.mdd2-se.mdd2, ymax=rg.mdd2+se.mdd2)) +
geom_linerange() +
geom_point() +
coord_flip(ylim=c(0, 1.25))
```

    ## Warning: Removed 4 rows containing missing values (geom_segment).

    ## Warning: Removed 4 rows containing missing values (geom_point).

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/metq_qc_ldsc_rg_mdd2-1.png)<!-- -->

Genetic covariance:

``` r
ggplot(meta_qc_ldsc, aes(x=reorder(paste(cohort, release), gencov.mdd2), y=gencov.mdd2, ymin=gencov.mdd2-gencov_se.mdd2, ymax=gencov.mdd2+gencov_se.mdd2)) +
geom_linerange() +
geom_point() +
coord_flip()
```

    ## Warning: Removed 3 rows containing missing values (geom_segment).

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/metq_qc_ldsc_gcov_mdd2-1.png)<!-- -->

# Heritabilities

Convert observed scale to liability scale with a range of population
prevalences:

``` r
meta_qc_h2_samp_prev <- 
meta_qc_ldsc %>%
transmute(cohort, release=str_sub(release, 2, -2), h2_obs, h2_obs_se) %>%
left_join(meta_qc_align %>%
            filter(ancestries == 'eur') %>%
            select(cohort, release, N_cases, N_controls),
          by=c('cohort', 'release')) %>%
mutate(samp_prev=N_cases / (N_cases + N_controls)) %>%
filter(!is.nan(h2_obs))

pop_prevs <- round(seq(0.05, 0.2, by=0.05), 2)

meta_qc_h2_liab <- 
plyr::adply(pop_prevs, 1, function(K) {
    meta_qc_h2_liab_pop <- 
    meta_qc_h2_samp_prev %>% mutate(K=K) %>%
    # normal distribution at pop prev threshold point
    mutate(zg=dnorm(qnorm(K))) %>%
    # call sample prevalence P
    mutate(P=samp_prev) %>%
    # liability scale h2
    mutate(h2_liab=h2_obs * K^2 * ( 1 - K)^2 / P / (1-P) / zg^2,
           h2_liab_se=h2_obs_se * K^2 * ( 1 - K)^2 / P / (1-P) / zg^2)
    return(meta_qc_h2_liab_pop)
}) %>%
mutate(h2_liab_low=h2_liab+qnorm(0.025)*h2_liab_se,
       h2_liab_upp=h2_liab+qnorm(0.975)*h2_liab_se)

ggplot(meta_qc_h2_liab %>% filter(K==0.1), aes(x=reorder(cohort, h2_liab), y=h2_liab, ymin=h2_liab_low, ymax=h2_liab_upp, group=cohort)) +
geom_hline(yintercept=c(0, 1), colour='gray') + 
geom_pointrange(position=position_dodge2()) +
coord_flip(ylim=c(-0.2, 1.5)) +
theme_minimal()
```

    ## Warning: Width not defined. Set with `position_dodge2(width = ?)`

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/meta_qc_h2_liab-1.png)<!-- -->

# Genetic covariance intercepts

Pairwise LDSC genetic covariance intercepts between all cohorts

``` r
meta_qc_ldsc_pairs <- read_tsv(snakemake@input$meta_qc_ldsc_pairs)
```

    ## 
    ## -- Column specification ----------------------------------------------------------------------------------------------------------
    ## cols(
    ##   cohort1 = col_character(),
    ##   subcohort1 = col_character(),
    ##   cohort2 = col_character(),
    ##   subcohort2 = col_character(),
    ##   ancestry1 = col_character(),
    ##   ancestry2 = col_character(),
    ##   rg = col_double(),
    ##   se = col_double(),
    ##   z = col_double(),
    ##   p = col_double(),
    ##   h2_obs = col_double(),
    ##   h2_obs_se = col_double(),
    ##   h2_int = col_double(),
    ##   h2_int_se = col_double(),
    ##   gcov_int = col_double(),
    ##   gcov_int_se = col_double()
    ## )

``` r
# calculate total sample sizes
meta_qc_samplesize <- meta_qc_align %>%
transmute(cohort, ancestries, release, N=N_cases+N_controls)

# sumstats appearing in MDD3 (exclude other
# sumstats like Wray 2018 that have had their 
# sumstats munged for other reasons)
meta_qc_ldsc_pairs_mdd3 <- meta_qc_ldsc_pairs  %>%
filter(!cohort1 %in% c('MDD29', 'PGC') & !cohort2 %in% c('MDD29', 'PGC')) %>%
left_join(meta_qc_samplesize %>% rename(cohortN1=N), by=c('cohort1'='cohort', 'subcohort1'='release', 'ancestry1'='ancestries')) %>%
left_join(meta_qc_samplesize %>% rename(cohortN2=N), by=c('cohort2'='cohort', 'subcohort2'='release', 'ancestry2'='ancestries'))
```

Histogram of intercepts

``` r
ggplot(meta_qc_ldsc_pairs_mdd3, aes(x=gcov_int)) +
geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 50 rows containing non-finite values (stat_bin).

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/meta_qc_ldcs_pairs_hist-1.png)<!-- -->

Pairwise intercepts

Show the largest and smallest intercepts. Exclude Wray2018 sumstats that
are also in the mix. Calculate distance from 0 in standard deviation
units to get an idea of the magnitude of departure from expectation
relative to the other intercepts.

``` r
meta_qc_ldsc_pairs_mdd3 %>%
select(cohort1, subcohort1, cohort2, subcohort2, gcov_int) %>%
mutate(SDs=abs(gcov_int)/sd(gcov_int, na.rm=T)) %>%
arrange(desc(gcov_int))
```

    ## # A tibble: 378 x 6
    ##    cohort1 subcohort1                 cohort2 subcohort2          gcov_int   SDs
    ##    <chr>   <chr>                      <chr>   <chr>                  <dbl> <dbl>
    ##  1 23andMe v7_2_202012                iPSYCH  2012_HRC              0.0219  3.46
    ##  2 MoBa    harvest12                  MoBa    rotterdam1            0.0213  3.36
    ##  3 MDD49   29w2_20w3_X28w2_19w3       tkda1   run1                  0.0212  3.34
    ##  4 HUNT    gp_hospital_metacarpa_201~ MoBa    harvest12             0.0203  3.20
    ##  5 HUNT    gp_hospital_metacarpa_201~ MVP     rel4icdDEP_Geno_20~   0.017   2.68
    ##  6 HUNT    gp_hospital_metacarpa_201~ iPSYCH  2015i_HRC             0.0165  2.60
    ##  7 UKBB    MD_glm_202107              deCODE  DEPALL_FINAL_WHEAD    0.0158  2.49
    ##  8 MDD49   29w2_20w3_X28w2_19w3       UKBB    MD_glm_202107         0.015   2.37
    ##  9 AGDS    202012                     iPSYCH  2012_HRC              0.0147  2.32
    ## 10 HUNT    gp_hospital_metacarpa_201~ iPSYCH  2012_HRC              0.0146  2.30
    ## # ... with 368 more rows

``` r
meta_qc_ldsc_pairs_mdd3 %>%
select(cohort1, subcohort1, cohort2, subcohort2, gcov_int) %>%
mutate(SDs=abs(gcov_int)/sd(gcov_int, na.rm=T)) %>%
arrange(gcov_int)
```

    ## # A tibble: 378 x 6
    ##    cohort1 subcohort1          cohort2 subcohort2          gcov_int   SDs
    ##    <chr>   <chr>               <chr>   <chr>                  <dbl> <dbl>
    ##  1 GERA    0915a_mds5          MoBa    harvest12            -0.0172  2.71
    ##  2 23andMe v7_2_202012         BASIC   202011               -0.0123  1.94
    ##  3 23andMe v7_2_202012         Airwave 0820                 -0.0121  1.91
    ##  4 SHARE   godartsshare_842021 iPSYCH  2012_HRC             -0.0116  1.83
    ##  5 EXCEED  202010              MoBa    harvest12            -0.0109  1.72
    ##  6 GERA    0915a_mds5          iPSYCH  2012_HRC             -0.0109  1.72
    ##  7 MoBa    harvest12           SHARE   godartsshare_842021  -0.0106  1.67
    ##  8 DBDS    FINAL202103         MoBa    rotterdam1           -0.0101  1.59
    ##  9 STAGE   MDDdx_saige         tkda1   run1                 -0.0089  1.40
    ## 10 GERA    0915a_mds5          PREFECT run1                 -0.0087  1.37
    ## # ... with 368 more rows

Estimate amount of sample overlap as
$N_S = g\_{\\mathrm{cov}\_\\mathrm{int}}\\sqrt{N_1N_2} / r\_\\mathrm{P}$,
where *r*<sub>P</sub> is the phenotypic correlation between the two MDD
phenotype (assume *r*<sub>P</sub> = 1)

``` r
rP <- 1
meta_qc_ldsc_pairs_mdd3_ns <-
meta_qc_ldsc_pairs_mdd3 %>%
mutate(Ns_factor=sqrt(cohortN1*cohortN2)/rP) %>%
mutate(Ns=gcov_int*Ns_factor,
       Ns_l95=(gcov_int+qnorm(0.025)*gcov_int_se)*Ns_factor,
       Ns_u95=(gcov_int+qnorm(0.975)*gcov_int_se)*Ns_factor,
       chisq=gcov_int^2/gcov_int_se^2) %>%
filter(!is.na(chisq)) %>%
mutate(qval=fdrtool(pchisq(chisq, df=1, lower.tail=F), statistic='pvalue', plot=FALSE)$qval) %>%
select(cohort1, subcohort1, cohort2, subcohort2, cohortN1, cohortN2, chisq, qval, Ns, Ns_l95, Ns_u95) %>%
mutate(Ns_pct=100*Ns/(cohortN1+cohortN2))
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

Sort (largest-to-smallest) by sample overlap

``` r
meta_qc_ldsc_pairs_mdd3_ns %>%
arrange(desc(Ns)) %>%
select(-chisq, -qval, -Ns_pct)
```

    ## # A tibble: 328 x 9
    ##    cohort1 subcohort1  cohort2 subcohort2 cohortN1 cohortN2    Ns  Ns_l95 Ns_u95
    ##    <chr>   <chr>       <chr>   <chr>         <dbl>    <dbl> <dbl>   <dbl>  <dbl>
    ##  1 23andMe v7_2_202012 UKBB    MD_glm_20~  1886830   361130 7512.  -5593. 20617.
    ##  2 23andMe v7_2_202012 MVP     rel4icdDE~  1886830   378614 6170. -12052. 24392.
    ##  3 23andMe v7_2_202012 iPSYCH  2012_HRC    1886830    41864 6155.   1032. 11278.
    ##  4 23andMe v7_2_202012 ESTBB   EstBB       1886830   126774 5429.  -1856. 12714.
    ##  5 23andMe v7_2_202012 BioVU   NoCov_SAI~  1886830    32480 2624.   -869.  6118.
    ##  6 ESTBB   EstBB       UKBB    MD_glm_20~   126774   361130 2568.    597.  4539.
    ##  7 23andMe v7_2_202012 GERA    0915a_mds5  1886830    45469 2548.  -1872.  6969.
    ##  8 MDD49   29w2_20w3_~ UKBB    MD_glm_20~    76180   361130 2488.    342.  4634.
    ##  9 HUNT    gp_hospita~ MVP     rel4icdDE~    54193   378614 2435.    442.  4428.
    ## 10 UKBB    MD_glm_202~ deCODE  DEPALL_FI~   361130    48000 2080.    738.  3422.
    ## # ... with 318 more rows

Sort (largest-to-smallest) by percentage of overlap to combined sample
size

``` r
meta_qc_ldsc_pairs_mdd3_ns %>%
arrange(desc(Ns_pct)) %>%
select(-cohortN1, -cohortN2, -chisq, -qval)
```

    ## # A tibble: 328 x 8
    ##    cohort1 subcohort1          cohort2 subcohort2        Ns Ns_l95 Ns_u95 Ns_pct
    ##    <chr>   <chr>               <chr>   <chr>          <dbl>  <dbl>  <dbl>  <dbl>
    ##  1 MoBa    harvest12           MoBa    rotterdam1     215.  102.     328.  1.06 
    ##  2 HUNT    gp_hospital_metaca~ iPSYCH  2015i_HRC      613.  220.    1006.  0.769
    ##  3 HUNT    gp_hospital_metaca~ MoBa    harvest12      491.  183.     800.  0.756
    ##  4 HUNT    gp_hospital_metaca~ iPSYCH  2012_HRC       695.   88.6   1302.  0.724
    ##  5 AGDS    202012              iPSYCH  2012_HRC       474.   56.9    891.  0.711
    ##  6 GenScot SCID_0721a          PBK     2020           173.   55.3    291.  0.611
    ##  7 STAGE   MDDdx_saige         lgic2   202011          92.4  16.2    168.  0.608
    ##  8 PREFECT run1                STAGE   MDDdx_saige     88.5   9.29   168.  0.605
    ##  9 iPSYCH  2012_HRC            iPSYCH  2015i_HRC      405.    1.70   808.  0.601
    ## 10 GERA    0915a_mds5          MDD49   29w2_20w3_X28~ 712.   54.6   1370.  0.585
    ## # ... with 318 more rows

Sort by *χ*<sup>2</sup> (Wald) test statistics

``` r
meta_qc_ldsc_pairs_mdd3_ns %>%
arrange(desc(chisq)) %>%
select(-cohortN1, -cohortN2)
```

    ## # A tibble: 328 x 10
    ##    cohort1 subcohort1 cohort2 subcohort2 chisq   qval    Ns Ns_l95 Ns_u95 Ns_pct
    ##    <chr>   <chr>      <chr>   <chr>      <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl>
    ##  1 MoBa    harvest12  MoBa    rotterdam1 14.0  0.0400  215.  102.    328.  1.06 
    ##  2 MDD49   29w2_20w3~ tkda1   run1       13.4  0.0400  228.  106.    350.  0.293
    ##  3 HUNT    gp_hospit~ MoBa    harvest12   9.75 0.134   491.  183.    800.  0.756
    ##  4 HUNT    gp_hospit~ iPSYCH  2015i_HRC   9.34 0.145   613.  220.   1006.  0.769
    ##  5 UKBB    MD_glm_20~ deCODE  DEPALL_FI~  9.23 0.148  2080.  738.   3422.  0.508
    ##  6 GenScot SCID_0721a PBK     2020        8.29 0.199   173.   55.3   291.  0.611
    ##  7 GenScot SCID_0721a iPSYCH  2012_HRC    7.46 0.246   237.   66.9   407.  0.489
    ##  8 GERA    0915a_mds5 MoBa    harvest12   7.45 0.246  -381. -655.   -108. -0.678
    ##  9 ESTBB   EstBB      UKBB    MD_glm_20~  6.52 0.328  2568.  597.   4539.  0.526
    ## 10 GERA    0915a_mds5 UKBB    MD_glm_20~  6.39 0.339  1717.  386.   3048.  0.422
    ## # ... with 318 more rows

Test for heterogeniety in covariance intercepts. Calculate *w*, the
inverse variance of each *g*<sub>covint</sub> standard error, then
calculate a weighted sum $\\hat{g\_\\mathrm{covint}}$

``` r
meta_qc_ldsc_pairs_mdd3 %>%
filter(!is.na(gcov_int)) %>%
mutate(w=1/gcov_int_se^2) %>%
mutate(gcov_int_hat=sum(w*gcov_int)/sum(w)) %>%
summarise(Q=sum(w*(gcov_int-gcov_int_hat)^2), k=n()) %>%
mutate(I2=(Q-(k-1))/Q)
```

    ## # A tibble: 1 x 3
    ##       Q     k    I2
    ##   <dbl> <int> <dbl>
    ## 1  382.   328 0.144

There is thus some heterogeneity in genetic covariance intercepts. We
also want to look at heterogeniety per-cohort. The LDSC intercepts were
calculated for each unique pair of cohorts. Expand this table to have
every ordering of each pair

``` r
# get unique cohort names from name and ancestry columns
meta_qc_ldsc_pairs_mdd3_named <-
meta_qc_ldsc_pairs_mdd3 %>%
mutate(cohort1=paste(cohort1, subcohort1, ancestry1, sep='.'),
       cohort2=paste(cohort2, subcohort2, ancestry2, sep='.')) %>%
mutate(pair_name=paste(cohort1, cohort2, sep=','))

# get all unique cohort names
cohorts <- unique(c(meta_qc_ldsc_pairs_mdd3_named$cohort1,
                    meta_qc_ldsc_pairs_mdd3_named$cohort2))

# ordered pairs of cohort we have data listings for
cohorts_pair_names <- unique(meta_qc_ldsc_pairs_mdd3_named$pair_name)

# create list of all pairs of cohorts
meta_qc_ldsc_pairs_mdd3_all <-
tibble(cohort1=cohorts, cohort2=cohorts) %>%
complete(cohort1, cohort2) %>%
filter(cohort1 != cohort2) %>%
mutate(pair_name12=paste(cohort1, cohort2, sep=','),
       pair_name21=paste(cohort2, cohort1, sep=',')) %>%
mutate(pair_name=case_when(pair_name12 %in% cohorts_pair_names ~ pair_name12,
                           pair_name21 %in% cohorts_pair_names ~ pair_name21,
                           TRUE ~ NA_character_)) %>%
select(cohort1, cohort2, pair_name) %>%
left_join(meta_qc_ldsc_pairs_mdd3_named %>% select(pair_name, gcov_int, gcov_int_se), by='pair_name')
```

Calculate *Q* and *I*<sup>2</sup> for each cohort

``` r
meta_qc_ldsc_pairs_mdd3_all_het <-
meta_qc_ldsc_pairs_mdd3_all %>%
group_by(cohort1) %>%
filter(!is.na(gcov_int)) %>%
mutate(w=1/gcov_int_se^2) %>%
mutate(gcov_int_hat=sum(w*gcov_int)/sum(w)) %>%
summarise(Q=sum(w*(gcov_int-gcov_int_hat)^2), k=n()) %>%
mutate(I2=(Q-(k-1))/Q) %>%
arrange(desc(I2))
meta_qc_ldsc_pairs_mdd3_all_het %>%
as.data.frame()
```

    ##                                    cohort1         Q  k           I2
    ## 1                       MoBa.harvest12.eur 46.485248 24  0.505219374
    ## 2                      GERA.0915a_mds5.eur 40.472200 24  0.431708681
    ## 3                           tkda1.run1.eur 25.959737 16  0.422182118
    ## 4                      iPSYCH.2012_HRC.eur 41.691342 27  0.376369315
    ## 5                      MoBa.rotterdam1.eur 36.826488 25  0.348295169
    ## 6                         BASIC.202011.eur 16.176164 12  0.319987114
    ## 7           MDD49.29w2_20w3_X28w2_19w3.eur 33.705239 25  0.287944534
    ## 8                         PREFECT.run1.eur 32.176958 24  0.285202788
    ## 9                    STAGE.MDDdx_saige.eur 32.221638 27  0.193088827
    ## 10 HUNT.gp_hospital_metacarpa_20190625.eur 29.350009 25  0.182283054
    ## 11                            PBK.2020.eur 24.857910 24  0.074741191
    ## 12                 23andMe.v7_2_202012.eur 26.693828 26  0.063453910
    ## 13                    iPSYCH.2015i_HRC.eur 25.534491 25  0.060094836
    ## 14                  GenScot.SCID_0721a.eur 24.407058 24  0.057649643
    ## 15                     ALSPAC.27022020.eur 21.032260 21  0.049079840
    ## 16                  UKBB.MD_glm_202107.eur 25.166148 25  0.046337967
    ## 17           deCODE.DEPALL_FINAL_WHEAD.eur 24.171284 25  0.007086266
    ## 18                        lgic2.202011.eur 23.675708 25 -0.013697258
    ## 19                         AGDS.202012.eur 23.000114 25 -0.043473109
    ## 20            BioVU.NoCov_SAIGE_051821.eur 20.979639 24 -0.096301062
    ## 21                         ESTBB.EstBB.eur 19.174763 25 -0.251645200
    ## 22           SHARE.godartsshare_842021.eur 17.333987 23 -0.269182878
    ## 23                    DBDS.FINAL202103.eur 19.449663 26 -0.285369320
    ## 24         MVP.rel4icdDEP_Geno_202109C.eur 18.591965 25 -0.290880236
    ## 25                       EXCEED.202010.eur 17.083763 25 -0.404842746
    ## 26                 FinnGen.R5_18032020.eur 16.545050 25 -0.450584902
    ## 27                        Airwave.0820.eur 16.252295 25 -0.476714548
    ## 28                      MoBa.harvest24.eur  3.921004  9 -1.040293511

``` r
# find range to fit all gcov_int values
gcov_int_min <- min(with(meta_qc_ldsc_pairs_mdd3, gcov_int-gcov_int_se), na.rm=TRUE)
gcov_int_max <- max(with(meta_qc_ldsc_pairs_mdd3, gcov_int+gcov_int_se), na.rm=TRUE)

ggplot(meta_qc_ldsc_pairs_mdd3_all %>% filter(cohort1 %in% meta_qc_ldsc_pairs_mdd3_all_het$cohort1[1:4]),
       aes(x=cohort2, y=gcov_int, ymin=gcov_int-gcov_int_se, ymax=gcov_int+gcov_int_se)) +
geom_point() +
geom_linerange() +
facet_grid(cols=vars(cohort1)) +
coord_flip(ylim=c(gcov_int_min-0.01, gcov_int_max+0.01))
```

    ## Warning: Removed 17 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_segment).

    ## Warning: Removed 3 rows containing missing values (geom_segment).

    ## Warning: Removed 11 rows containing missing values (geom_segment).

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/plot_cohort_gcov_int-1.png)<!-- -->

``` r
ggplot(meta_qc_ldsc_pairs_mdd3_all %>% filter(cohort1 %in% meta_qc_ldsc_pairs_mdd3_all_het$cohort1[5:8]),
   aes(x=cohort2, y=gcov_int, ymin=gcov_int-gcov_int_se, ymax=gcov_int+gcov_int_se)) +
geom_point() +
geom_linerange() +
facet_grid(cols=vars(cohort1)) +
coord_flip(ylim=c(gcov_int_min-0.01, gcov_int_max+0.01))
```

    ## Warning: Removed 22 rows containing missing values (geom_point).

    ## Warning: Removed 15 rows containing missing values (geom_segment).

    ## Warning: Removed 2 rows containing missing values (geom_segment).

    ## Warning: Removed 2 rows containing missing values (geom_segment).

    ## Warning: Removed 3 rows containing missing values (geom_segment).

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/plot_cohort_gcov_int2-1.png)<!-- -->

## Clustering

Cluster based on similarity in genetic covariance intercepts.

``` r
# get a list of all cohort/subcohort names

subcohorts <- 
bind_rows(
select(meta_qc_ldsc_pairs_mdd3, cohort=cohort1, subcohort=subcohort1),
select(meta_qc_ldsc_pairs_mdd3, cohort=cohort2, subcohort=subcohort2)
) %>%
distinct()
```

Make a matrix:

``` r
cohort_names <- 
subcohorts %>%
transmute(cohort=paste(cohort, subcohort, sep='.')) %>%
pull(cohort)

gcov_int_mat <- diag(length(cohort_names))
dimnames(gcov_int_mat) <- list(cohort_names, cohort_names)

for(i in seq.int(nrow(meta_qc_ldsc_pairs_mdd3))) {
    cohort1 <- meta_qc_ldsc_pairs_mdd3$cohort1[i]
    cohort2 <- meta_qc_ldsc_pairs_mdd3$cohort2[i]
    subcohort1 <- meta_qc_ldsc_pairs_mdd3$subcohort1[i]
    subcohort2 <- meta_qc_ldsc_pairs_mdd3$subcohort2[i]
    cohort_name1 <- paste(cohort1, subcohort1, sep='.')
    cohort_name2 <- paste(cohort2, subcohort2, sep='.')

    # if gcov_int is NA, substitute 0 (average)

    gcov_int <- coalesce(meta_qc_ldsc_pairs_mdd3$gcov_int[i], 0)

    gcov_int_mat[cohort_name1,cohort_name2] <- gcov_int 
    gcov_int_mat[cohort_name2,cohort_name1] <- gcov_int
}
```

Convert intercepts into a dissimilarity matrix. For intercepts, `1 ==`
same, `0 ==` average, `<0 ==` different. For dissimilarity, `0 ==` same
and larger values `==` more dissimilar.

``` r
plot(hclust(as.dist(1-gcov_int_mat)))
```

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/gcov_int_dist-1.png)<!-- -->

``` r
corrplot(gcov_int_mat, is.corr=FALSE, diag=FALSE)
```

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/gcov_int_mat_corr-1.png)<!-- -->

## Genetic correlations

``` r
rg_mat <- diag(length(cohort_names))
dimnames(rg_mat) <- list(cohort_names, cohort_names)
rg_se_mat <- diag(x=0, nrow=length(cohort_names), ncol=length(cohort_names))
dimnames(rg_se_mat) <- list(cohort_names, cohort_names)

for(i in seq.int(nrow(meta_qc_ldsc_pairs_mdd3))) {
    cohort1 <- meta_qc_ldsc_pairs_mdd3$cohort1[i]
    cohort2 <- meta_qc_ldsc_pairs_mdd3$cohort2[i]
    subcohort1 <- meta_qc_ldsc_pairs_mdd3$subcohort1[i]
    subcohort2 <- meta_qc_ldsc_pairs_mdd3$subcohort2[i]
    cohort_name1 <- paste(cohort1, subcohort1, sep='.')
    cohort_name2 <- paste(cohort2, subcohort2, sep='.')

    # if rg is NA, substitute 0 (average)

    rg <- meta_qc_ldsc_pairs_mdd3$rg[i]
    rg_se <- meta_qc_ldsc_pairs_mdd3$se[i]

    rg_mat[cohort_name1,cohort_name2] <- rg 
    rg_mat[cohort_name2,cohort_name1] <- rg
    rg_se_mat[cohort_name1,cohort_name2] <- rg_se
    rg_se_mat[cohort_name2,cohort_name1] <- rg_se
}
```

``` r
rg_mat_1 <- rg_mat
rg_mat_1[which(rg_mat_1 > 1)] <- 1
rg_mat_1[which(rg_mat_1 < -1)] <- -1

rg_lowCI <- rg_mat + qnorm(0.025)*rg_se_mat
rg_uppCI <- rg_mat + qnorm(0.975)*rg_se_mat

rg_lowCI_1 <- rg_lowCI
rg_lowCI_1[which(rg_lowCI_1 > 1)] <- 1
rg_lowCI_1[which(rg_lowCI_1 < -1)] <- -1

rg_uppCI_1 <- rg_uppCI
rg_uppCI_1[which(rg_uppCI_1 > 1)] <- 1
rg_uppCI_1[which(rg_uppCI_1 < -1)] <- -1




has_rg_idx <- rowSums(!is.na(rg_mat_1)) > 1

corrplot.mixed(rg_mat_1[has_rg_idx,has_rg_idx],
               lowCI=rg_lowCI_1[has_rg_idx,has_rg_idx],
               uppCI=rg_uppCI_1[has_rg_idx,has_rg_idx],
               na.label='.', tl.pos='lt', diag='u', plotCI='rect', number.cex=0.75)
```

![](/home/madams/projects/mdd-meta/docs/metaqc_files/figure-gfm/rg_mat_corr-1.png)<!-- -->
