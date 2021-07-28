Structure Meta-Analysis with GenomicSEM
================

Summary statistics were grouped based on phenotype, then meta-analysed:

-   `clin`: Clinical assessment
-   `ehr`: Electronic health records
-   `quest`: Questionnaire
-   `self`: Single-item self report

The LDSC covariance structure of the four MDD phenotype was calcuted
using [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM).

``` r
library(GenomicSEM)
library(readr)
library(corrplot)
library(dplyr)
```

Read in the covariance structure and the LDSC

``` r
mdd_covstruct <- dget(snakemake@input$covstruct)
```

LDSC statistics

``` r
ldsc_table <- read_tsv(snakemake@input$ldsc_table)
```

    ## Rows: 4 Columns: 14

    ## ── Column specification ─────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (2): pheno, ancestries
    ## dbl (12): N_cases, N_controls, sample_prev, pop_prev, LambdaGC, MeanChiSq, L...

    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
knitr::kable(ldsc_table)
```

| pheno | ancestries | N\_cases | N\_controls | sample\_prev | pop\_prev | LambdaGC | MeanChiSq | LambdaGCldsc | ldsc\_intercept | h2\_obs | h2\_se\_obs | h2\_liab | h2\_se\_liab |
|:------|:-----------|---------:|------------:|-------------:|----------:|---------:|----------:|-------------:|----------------:|--------:|------------:|---------:|-------------:|
| clin  | eur        |    30722 |       59954 |      0.33880 |      0.15 |    1.124 |    1.1484 |       1.1490 |          1.0214 |  0.0680 |      0.0071 |   0.0908 |       0.0095 |
| ehr   | eur        |   310522 |      876421 |      0.26160 |      0.15 |    1.444 |    1.7548 |       1.5807 |          1.0279 |  0.0313 |      0.0011 |   0.0485 |       0.0018 |
| quest | eur        |    55519 |      346777 |      0.13800 |      0.15 |    1.194 |    1.2725 |       1.2266 |          1.0156 |  0.0328 |      0.0018 |   0.0824 |       0.0045 |
| self  | eur        |   114992 |     1789651 |      0.06037 |      0.15 |    1.453 |    1.9058 |       1.6485 |          0.9803 |  0.0240 |      0.0009 |   0.1265 |       0.0049 |

Genetic correlations

``` r
corrplot(cov2cor(mdd_covstruct$S), method='number')
```

![](gsem_files/figure-gfm/rg-1.png)<!-- -->

## Common factor model

``` r
common.model <- "A =~ NA*clin + ehr + quest + self
A ~~ 1*A"

common.fit <- usermodel(covstruc=mdd_covstruct, estimation='DWLS', model=common.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.182

``` r
knitr::kable(common.fit$modelfit)
```

|     |    chisq |  df |  p\_chisq |      AIC | CFI |      SRMR |
|:----|---------:|----:|----------:|---------:|----:|----------:|
| df  | 0.612646 |   2 | 0.7361488 | 16.61265 |   1 | 0.0075507 |

``` r
common.fit.results <- common.fit$results %>%
    mutate(Unstand_SE=as.numeric(Unstand_SE),
           STD_Genotype_SE=as.numeric(STD_Genotype_SE),
           p_value=if_else(p_value == '< 5e-300', true=0, false=as.numeric(p_value)))
```

    ## Warning in replace_with(out, !condition, false, fmt_args(~false), glue("length
    ## of {fmt_args(~condition)}")): NAs introduced by coercion

``` r
knitr::kable(common.fit.results, digits=4)
```

|     | lhs   | op   | rhs   | Unstand\_Est | Unstand\_SE | STD\_Genotype | STD\_Genotype\_SE | STD\_All | p\_value |
|:----|:------|:-----|:------|-------------:|------------:|--------------:|------------------:|---------:|---------:|
| 1   | A     | =\~  | clin  |       0.2930 |      0.0121 |        0.9724 |            0.0403 |   0.9724 |   0.0000 |
| 2   | A     | =\~  | ehr   |       0.2074 |      0.0044 |        0.9417 |            0.0200 |   0.9417 |   0.0000 |
| 3   | A     | =\~  | quest |       0.2758 |      0.0082 |        0.9609 |            0.0285 |   0.9609 |   0.0000 |
| 4   | A     | =\~  | self  |       0.2945 |      0.0087 |        0.8282 |            0.0245 |   0.8282 |   0.0000 |
| 6   | clin  | \~\~ | clin  |       0.0049 |      0.0078 |        0.0544 |            0.0854 |   0.0544 |   0.5246 |
| 7   | ehr   | \~\~ | ehr   |       0.0055 |      0.0015 |        0.1132 |            0.0318 |   0.1132 |   0.0004 |
| 8   | quest | \~\~ | quest |       0.0063 |      0.0043 |        0.0767 |            0.0520 |   0.0767 |   0.1400 |
| 9   | self  | \~\~ | self  |       0.0397 |      0.0037 |        0.3140 |            0.0293 |   0.3140 |   0.0000 |
| 5   | A     | \~\~ | A     |       1.0000 |          NA |        1.0000 |                NA |   1.0000 |       NA |

## GWAS-by-subtraction

Fit a [GWAS-by-subtraction](https://rpubs.com/MichelNivard/565885) to
decompose variance into that of minimally-phenotyped depression shared
with maximal phenotypes versus specific to maximal phenotypes.

``` r
md.model <- "MD =~ NA*self + clin + ehr + quest
MDD =~ NA*clin + ehr + quest

MD ~~ 1*MD
MDD ~~ 1*MDD
MD ~~ 0*MDD

self ~~ 0*clin + 0*ehr + 0*quest
self ~~ 0*self"
```

``` r
md.fit <- usermodel(covstruc=mdd_covstruct, estimation='DWLS', model=md.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##    0.18 
    ## [1] "Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)"

``` r
md.fit.results <- md.fit$results %>%
mutate(Unstand_SE=as.numeric(Unstand_SE),
       STD_Genotype_SE=as.numeric(STD_Genotype_SE),
       p_value=if_else(p_value == '< 5e-300', true=0, false=as.numeric(p_value)))
```

    ## Warning in replace_with(out, !condition, false, fmt_args(~false), glue("length
    ## of {fmt_args(~condition)}")): NAs introduced by coercion

``` r
knitr::kable(md.fit.results, digits=4)
```

|     | lhs   | op   | rhs   | Unstand\_Est | Unstand\_SE | STD\_Genotype | STD\_Genotype\_SE | STD\_All | p\_value |
|:----|:------|:-----|:------|-------------:|------------:|--------------:|------------------:|---------:|---------:|
| 6   | MD    | =\~  | self  |       0.3556 |      0.0069 |        1.0000 |            0.0194 |   1.0000 |   0.0000 |
| 3   | MD    | =\~  | clin  |       0.2459 |      0.0119 |        0.8162 |            0.0393 |   0.8162 |   0.0000 |
| 4   | MD    | =\~  | ehr   |       0.1707 |      0.0045 |        0.7752 |            0.0203 |   0.7752 |   0.0000 |
| 5   | MD    | =\~  | quest |       0.2286 |      0.0082 |        0.7966 |            0.0287 |   0.7966 |   0.0000 |
| 9   | MDD   | =\~  | clin  |       0.1523 |      0.0208 |        0.5056 |            0.0691 |   0.5056 |   0.0000 |
| 10  | MDD   | =\~  | ehr   |       0.1238 |      0.0135 |        0.5622 |            0.0614 |   0.5622 |   0.0000 |
| 11  | MDD   | =\~  | quest |       0.1502 |      0.0162 |        0.5235 |            0.0563 |   0.5235 |   0.0000 |
| 1   | clin  | \~\~ | clin  |       0.0071 |      0.0082 |        0.0782 |            0.0899 |   0.0782 |   0.3843 |
| 2   | ehr   | \~\~ | ehr   |       0.0040 |      0.0028 |        0.0831 |            0.0582 |   0.0831 |   0.1534 |
| 13  | quest | \~\~ | quest |       0.0075 |      0.0052 |        0.0913 |            0.0635 |   0.0913 |   0.1506 |
| 7   | MD    | \~\~ | MD    |       1.0000 |          NA |        1.0000 |                NA |   1.0000 |       NA |
| 12  | MDD   | \~\~ | MDD   |       1.0000 |          NA |        1.0000 |                NA |   1.0000 |       NA |
| 8   | MD    | \~\~ | MDD   |       0.0000 |          NA |        0.0000 |                NA |   0.0000 |       NA |
| 14  | self  | \~\~ | clin  |       0.0000 |          NA |        0.0000 |                NA |   0.0000 |       NA |
| 15  | self  | \~\~ | ehr   |       0.0000 |          NA |        0.0000 |                NA |   0.0000 |       NA |
| 16  | self  | \~\~ | quest |       0.0000 |          NA |        0.0000 |                NA |   0.0000 |       NA |
| 17  | self  | \~\~ | self  |       0.0000 |          NA |        0.0000 |                NA |   0.0000 |       NA |
