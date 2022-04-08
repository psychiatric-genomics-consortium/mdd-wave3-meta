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
library(miamiplot)
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
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (2): pheno, ancestries
    ## dbl (12): N_cases, N_controls, sample_prev, pop_prev, LambdaGC, MeanChiSq, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
knitr::kable(ldsc_table)
```

| pheno   | ancestries | N_cases | N_controls | sample_prev | pop_prev | LambdaGC | MeanChiSq | LambdaGCldsc | ldsc_intercept | h2_obs | h2_se_obs | h2_liab | h2_se_liab |
|:--------|:-----------|--------:|-----------:|------------:|---------:|---------:|----------:|-------------:|---------------:|-------:|----------:|--------:|-----------:|
| Clin    | eur        |   26778 |      51857 |         0.5 |     0.15 |    1.104 |    1.1353 |       1.1322 |         1.0126 | 0.0956 |    0.0099 |  0.1143 |     0.0118 |
| EHR     | eur        |  311831 |    1180439 |         0.5 |     0.15 |    1.483 |    1.8260 |       1.6305 |         1.0323 | 0.0490 |    0.0018 |  0.0586 |     0.0021 |
| Quest   | eur        |   70887 |     339597 |         0.5 |     0.15 |    1.219 |    1.3219 |       1.2716 |         1.0167 | 0.0705 |    0.0036 |  0.0843 |     0.0043 |
| SelfRep | eur        |  114992 |    1789651 |         0.5 |     0.15 |    1.580 |    1.9628 |       1.7091 |         1.0034 | 0.1080 |    0.0041 |  0.1292 |     0.0050 |

Genetic correlations

``` r
corrplot(cov2cor(mdd_covstruct$S), method = "number")
```

![](/gpfs/igmmfs01/eddie/GenScotDepression/madams23/projects/mdd-meta/docs/gsem_files/figure-gfm/rg-1.png)<!-- -->

## Common factor model

``` r
common.model <- "A =~ NA*Clin + EHR + Quest + SelfRep
A ~~ 1*A"

common.fit <- usermodel(covstruc = mdd_covstruct,
                        estimation = "DWLS",
                        model = common.model,
                        imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.745

``` r
knitr::kable(common.fit$modelfit)
```

|     |    chisq |  df |   p_chisq |      AIC |       CFI |     SRMR |
|:----|---------:|----:|----------:|---------:|----------:|---------:|
| df  | 2.100548 |   2 | 0.3498419 | 18.10055 | 0.9999705 | 0.011537 |

``` r
common.fit.results <- common.fit$results %>%
    mutate(Unstand_SE = as.numeric(Unstand_SE),
           STD_Genotype_SE = as.numeric(STD_Genotype_SE),
           p_value=if_else(p_value == "< 5e-300",
                           true = 0,
                           false = as.numeric(p_value)))
```

    ## Warning in replace_with(out, !condition, false, "`false`", "length of
    ## `condition`"): NAs introduced by coercion

``` r
knitr::kable(select(common.fit.results,
                    lhs, op, rhs,
                    STD_Genotype, STD_Genotype_SE,
                    p_value),
                    digits = 4)
```

|     | lhs     | op   | rhs     | STD_Genotype | STD_Genotype_SE | p_value |
|:----|:--------|:-----|:--------|-------------:|----------------:|--------:|
| 1   | A       | =\~  | Clin    |       0.9178 |          0.0349 |  0.0000 |
| 2   | A       | =\~  | EHR     |       0.9208 |          0.0212 |  0.0000 |
| 3   | A       | =\~  | Quest   |       0.9532 |          0.0285 |  0.0000 |
| 4   | A       | =\~  | SelfRep |       0.8454 |          0.0295 |  0.0000 |
| 6   | Clin    | \~\~ | Clin    |       0.1577 |          0.0912 |  0.0839 |
| 7   | EHR     | \~\~ | EHR     |       0.1522 |          0.0320 |  0.0000 |
| 8   | Quest   | \~\~ | Quest   |       0.0914 |          0.0523 |  0.0804 |
| 9   | SelfRep | \~\~ | SelfRep |       0.2853 |          0.0425 |  0.0000 |
| 5   | A       | \~\~ | A       |       1.0000 |              NA |      NA |

Model-implied correlations

``` r
cov2cor(common.fit$resid_cov$Model)
```

    ##              Clin       EHR     Quest  SelfRep
    ## Clin    1.0000000 0.8450678 0.8748459 0.775885
    ## EHR     0.8450678 1.0000000 0.8776829 0.778401
    ## Quest   0.8748459 0.8776829 1.0000000 0.805830
    ## SelfRep 0.7758850 0.7784010 0.8058300 1.000000

### Common factor GWAS

Miami plot of P values and Q values

``` r
common.gwas <- read_tsv(snakemake@input$commonfactor_eur)
```

    ## Rows: 6399123 Columns: 19
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (6): SNP, A1, A2, lhs, op, rhs
    ## dbl (13): CHR, BP, MAF, i, est, se_c, Z_Estimate, Pval_Estimate, Q, Q_df, Q_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
p.gwas <- common.gwas %>% transmute(rsid=SNP, chr=CHR, pos=BP, pval=Pval_Estimate, P="P")
q.gwas <- common.gwas %>% transmute(rsid=SNP, chr=CHR, pos=BP, pval=Q_pval, P="Q")

ggmiami(data=bind_rows(p.gwas, q.gwas) %>% filter(pval <= 0.001) %>% as.data.frame(),
        split_by="P", split_at="P", p="pval",
        upper_ylab="Common Factor GWAS",
        lower_ylab="Heterogeneity Q")
```

![](/gpfs/igmmfs01/eddie/GenScotDepression/madams23/projects/mdd-meta/docs/gsem_files/figure-gfm/gwas-1.png)<!-- -->

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
    ##   0.756 
    ## [1] "Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)"

``` r
md.fit.results <- md.fit$results %>%
mutate(Unstand_SE=as.numeric(Unstand_SE),
       STD_Genotype_SE=as.numeric(STD_Genotype_SE),
       p_value=if_else(p_value == '< 5e-300', true=0, false=as.numeric(p_value)))
```

    ## Warning in replace_with(out, !condition, false, "`false`", "length of
    ## `condition`"): NAs introduced by coercion

``` r
knitr::kable(md.fit.results, digits=4)
```

|     | lhs     | op   | rhs     | Unstand_Est | Unstand_SE | STD_Genotype | STD_Genotype_SE | STD_All | p_value |
|:----|:--------|:-----|:--------|------------:|-----------:|-------------:|----------------:|--------:|--------:|
| 10  | MIN     | =\~  | SelfRep |      0.3595 |     0.0069 |       1.0000 |          0.0192 |  1.0000 |  0.0000 |
| 7   | MIN     | =\~  | Clin    |      0.2622 |     0.0133 |       0.7753 |          0.0393 |  0.7753 |  0.0000 |
| 8   | MIN     | =\~  | EHR     |      0.1864 |     0.0056 |       0.7700 |          0.0231 |  0.7700 |  0.0000 |
| 9   | MIN     | =\~  | Quest   |      0.2377 |     0.0082 |       0.8186 |          0.0281 |  0.8186 |  0.0000 |
| 3   | MAX     | =\~  | Clin    |      0.1632 |     0.0305 |       0.4827 |          0.0903 |  0.4827 |  0.0000 |
| 4   | MAX     | =\~  | EHR     |      0.1347 |     0.0167 |       0.5566 |          0.0691 |  0.5566 |  0.0000 |
| 5   | MAX     | =\~  | Quest   |      0.1289 |     0.0208 |       0.4438 |          0.0718 |  0.4438 |  0.0000 |
| 1   | Clin    | \~\~ | Clin    |      0.0190 |     0.0117 |       0.1659 |          0.1019 |  0.1659 |  0.1037 |
| 2   | EHR     | \~\~ | EHR     |      0.0057 |     0.0038 |       0.0974 |          0.0648 |  0.0974 |  0.1331 |
| 13  | Quest   | \~\~ | Quest   |      0.0112 |     0.0047 |       0.1330 |          0.0557 |  0.1330 |  0.0170 |
| 12  | MIN     | \~\~ | MIN     |      1.0000 |         NA |       1.0000 |              NA |  1.0000 |      NA |
| 6   | MAX     | \~\~ | MAX     |      1.0000 |         NA |       1.0000 |              NA |  1.0000 |      NA |
| 11  | MIN     | \~\~ | MAX     |      0.0000 |         NA |       0.0000 |              NA |  0.0000 |      NA |
| 14  | SelfRep | \~\~ | Clin    |      0.0000 |         NA |       0.0000 |              NA |  0.0000 |      NA |
| 15  | SelfRep | \~\~ | EHR     |      0.0000 |         NA |       0.0000 |              NA |  0.0000 |      NA |
| 16  | SelfRep | \~\~ | Quest   |      0.0000 |         NA |       0.0000 |              NA |  0.0000 |      NA |
| 17  | SelfRep | \~\~ | SelfRep |      0.0000 |         NA |       0.0000 |              NA |  0.0000 |      NA |
