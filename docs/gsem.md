Structure Meta-Analysis with GenomicSEM
================

Summary statistics were grouped based on phenotype, then meta-analysed:

-   `Clin`: Clinical assessment
-   `EHR`: Electronic health records
-   `Quest`: Questionnaire
-   `SelfRep`: Single-item self report

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

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   pheno = col_character(),
    ##   ancestries = col_character(),
    ##   N_cases = col_double(),
    ##   N_controls = col_double(),
    ##   sample_prev = col_double(),
    ##   pop_prev = col_double(),
    ##   LambdaGC = col_double(),
    ##   MeanChiSq = col_double(),
    ##   LambdaGCldsc = col_double(),
    ##   ldsc_intercept = col_double(),
    ##   h2_obs = col_double(),
    ##   h2_se_obs = col_double(),
    ##   h2_liab = col_double(),
    ##   h2_se_liab = col_double()
    ## )

``` r
knitr::kable(ldsc_table)
```

| pheno   | ancestries | N_cases | N_controls | sample_prev | pop_prev | LambdaGC | MeanChiSq | LambdaGCldsc | ldsc_intercept | h2_obs | h2_se_obs | h2_liab | h2_se_liab |
|:--------|:-----------|--------:|-----------:|------------:|---------:|---------:|----------:|-------------:|---------------:|-------:|----------:|--------:|-----------:|
| Clin    | eur        |   21654 |      44593 |         0.5 |     0.15 |    1.097 |    1.1087 |       1.1116 |         1.0227 | 0.0811 |    0.0110 |  0.0970 |     0.0131 |
| EHR     | eur        |  311491 |     877110 |         0.5 |     0.15 |    1.516 |    1.7737 |       1.5964 |         1.0363 | 0.0532 |    0.0019 |  0.0637 |     0.0023 |
| Quest   | eur        |   70887 |     339597 |         0.5 |     0.15 |    1.230 |    1.3224 |       1.2716 |         1.0161 | 0.0708 |    0.0036 |  0.0847 |     0.0043 |
| SelfRep | eur        |  114992 |    1789651 |         0.5 |     0.15 |    1.645 |    1.9633 |       1.7098 |         1.0042 | 0.1080 |    0.0042 |  0.1291 |     0.0050 |

Genetic correlations

``` r
corrplot(cov2cor(mdd_covstruct$S), method='number')
```

![](/gsem_files/figure-gfm/rg-1.png)<!-- -->

## Common factor model

``` r
common.model <- "A =~ NA*Clin + EHR + Quest + SelfRep
A ~~ 1*A"

common.fit <- usermodel(covstruc=mdd_covstruct, estimation='DWLS', model=common.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.567

``` r
knitr::kable(common.fit$modelfit)
```

|     |    chisq |  df |   p_chisq |     AIC | CFI |      SRMR |
|:----|---------:|----:|----------:|--------:|----:|----------:|
| df  | 1.783197 |   2 | 0.4099997 | 17.7832 |   1 | 0.0142189 |

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

|     | lhs     | op   | rhs     | Unstand_Est | Unstand_SE | STD_Genotype | STD_Genotype_SE | STD_All | p_value |
|:----|:--------|:-----|:--------|------------:|-----------:|-------------:|----------------:|--------:|--------:|
| 1   | A       | =\~  | Clin    |      0.3101 |     0.0121 |       0.9957 |          0.0388 |  0.9957 |  0.0000 |
| 2   | A       | =\~  | EHR     |      0.2360 |     0.0052 |       0.9353 |          0.0207 |  0.9353 |  0.0000 |
| 3   | A       | =\~  | Quest   |      0.2751 |     0.0083 |       0.9452 |          0.0286 |  0.9452 |  0.0000 |
| 4   | A       | =\~  | SelfRep |      0.3075 |     0.0101 |       0.8556 |          0.0281 |  0.8556 |  0.0000 |
| 6   | Clin    | \~\~ | Clin    |      0.0008 |     0.0114 |       0.0085 |          0.1173 |  0.0085 |  0.9420 |
| 7   | EHR     | \~\~ | EHR     |      0.0080 |     0.0021 |       0.1252 |          0.0329 |  0.1252 |  0.0001 |
| 8   | Quest   | \~\~ | Quest   |      0.0090 |     0.0047 |       0.1065 |          0.0557 |  0.1065 |  0.0560 |
| 9   | SelfRep | \~\~ | SelfRep |      0.0346 |     0.0055 |       0.2680 |          0.0426 |  0.2680 |  0.0000 |
| 5   | A       | \~\~ | A       |      1.0000 |         NA |       1.0000 |              NA |  1.0000 |      NA |

## GWAS-by-subtraction

Fit a [GWAS-by-subtraction](https://rpubs.com/MichelNivard/565885) to
decompose variance into that of minimally-phenotyped depression shared
with maximal phenotypes versus specific to maximal phenotypes.

``` r
md.model <- "MIN =~ NA*SelfRep + Clin + EHR + Quest
MAX =~ NA*Clin + EHR + Quest

MIN ~~ 1*MIN
MAX ~~ 1*MAX
MIN ~~ 0*MAX

SelfRep ~~ 0*Clin + 0*EHR + 0*Quest
SelfRep ~~ 0*SelfRep"
```

``` r
md.fit <- usermodel(covstruc=mdd_covstruct, estimation='DWLS', model=md.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.326 
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

|     | lhs     | op   | rhs     | Unstand_Est | Unstand_SE | STD_Genotype | STD_Genotype_SE | STD_All | p_value |
|:----|:--------|:-----|:--------|------------:|-----------:|-------------:|----------------:|--------:|--------:|
| 10  | MIN     | =\~  | SelfRep |      0.3594 |     0.0069 |       1.0000 |          0.0192 |  1.0000 |  0.0000 |
| 7   | MIN     | =\~  | Clin    |      0.2674 |     0.0139 |       0.8586 |          0.0446 |  0.8586 |  0.0000 |
| 8   | MIN     | =\~  | EHR     |      0.1995 |     0.0060 |       0.7909 |          0.0237 |  0.7909 |  0.0000 |
| 9   | MIN     | =\~  | Quest   |      0.2385 |     0.0081 |       0.8193 |          0.0279 |  0.8193 |  0.0000 |
| 3   | MAX     | =\~  | Clin    |      0.1490 |     0.0315 |       0.4785 |          0.1011 |  0.4785 |  0.0000 |
| 4   | MAX     | =\~  | EHR     |      0.1439 |     0.0220 |       0.5705 |          0.0871 |  0.5705 |  0.0000 |
| 5   | MAX     | =\~  | Quest   |      0.1224 |     0.0231 |       0.4205 |          0.0793 |  0.4205 |  0.0000 |
| 1   | Clin    | \~\~ | Clin    |      0.0033 |     0.0120 |       0.0339 |          0.1241 |  0.0339 |  0.7847 |
| 2   | EHR     | \~\~ | EHR     |      0.0031 |     0.0054 |       0.0490 |          0.0847 |  0.0490 |  0.5626 |
| 13  | Quest   | \~\~ | Quest   |      0.0129 |     0.0052 |       0.1519 |          0.0617 |  0.1519 |  0.0139 |
| 12  | MIN     | \~\~ | MIN     |      1.0000 |         NA |       1.0000 |              NA |  1.0000 |      NA |
| 6   | MAX     | \~\~ | MAX     |      1.0000 |         NA |       1.0000 |              NA |  1.0000 |      NA |
| 11  | MIN     | \~\~ | MAX     |      0.0000 |         NA |       0.0000 |              NA |  0.0000 |      NA |
| 14  | SelfRep | \~\~ | Clin    |      0.0000 |         NA |       0.0000 |              NA |  0.0000 |      NA |
| 15  | SelfRep | \~\~ | EHR     |      0.0000 |         NA |       0.0000 |              NA |  0.0000 |      NA |
| 16  | SelfRep | \~\~ | Quest   |      0.0000 |         NA |       0.0000 |              NA |  0.0000 |      NA |
| 17  | SelfRep | \~\~ | SelfRep |      0.0000 |         NA |       0.0000 |              NA |  0.0000 |      NA |
