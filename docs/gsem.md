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
```

    ## Warning: replacing previous import 'gdata::nobs' by 'lavaan::nobs' when loading
    ## 'GenomicSEM'

    ## Warning: replacing previous import 'gdata::last' by 'data.table::last' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'gdata::first' by 'data.table::first' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'gdata::env' by 'R.utils::env' when loading
    ## 'GenomicSEM'

    ## Warning: replacing previous import 'gdata::resample' by 'R.utils::resample' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'data.table::last' by 'dplyr::last' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'plyr::summarize' by 'dplyr::summarize' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'plyr::mutate' by 'dplyr::mutate' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'plyr::id' by 'dplyr::id' when loading
    ## 'GenomicSEM'

    ## Warning: replacing previous import 'plyr::arrange' by 'dplyr::arrange' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'plyr::summarise' by 'dplyr::summarise' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'data.table::first' by 'dplyr::first' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'plyr::rename' by 'dplyr::rename' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'plyr::desc' by 'dplyr::desc' when loading
    ## 'GenomicSEM'

    ## Warning: replacing previous import 'plyr::count' by 'dplyr::count' when loading
    ## 'GenomicSEM'

    ## Warning: replacing previous import 'gdata::combine' by 'dplyr::combine' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'data.table::between' by 'dplyr::between'
    ## when loading 'GenomicSEM'

    ## Warning: replacing previous import 'plyr::failwith' by 'dplyr::failwith' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'Rcpp::prompt' by 'utils::prompt' when
    ## loading 'GenomicSEM'

    ## Warning: replacing previous import 'Rcpp::.DollarNames' by 'utils::.DollarNames'
    ## when loading 'GenomicSEM'

    ## Warning: replacing previous import 'Matrix::tail' by 'utils::tail' when loading
    ## 'GenomicSEM'

    ## Warning: replacing previous import 'gdata::object.size' by 'utils::object.size'
    ## when loading 'GenomicSEM'

    ## Warning: replacing previous import 'R.utils::timestamp' by 'utils::timestamp'
    ## when loading 'GenomicSEM'

    ## Warning: replacing previous import 'Matrix::head' by 'utils::head' when loading
    ## 'GenomicSEM'

``` r
library(readr)
library(corrplot)
```

    ## corrplot 0.90 loaded

Read in the covariance structure and the LDSC

``` r
mdd_covstruct <- dget(snakemake@input$covstruct)
```

LDSC statistics

``` r
ldsc_table <- read_tsv(snakemake@input$ldsc_table)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────────
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

| pheno | ancestries | N_cases | N_controls | sample_prev | pop_prev | LambdaGC | MeanChiSq | LambdaGCldsc | ldsc_intercept | h2_obs | h2_se_obs | h2_liab | h2_se_liab |
|:------|:-----------|--------:|-----------:|------------:|---------:|---------:|----------:|-------------:|---------------:|-------:|----------:|--------:|-----------:|
| clin  | eur        |   30722 |      59954 |     0.33880 |     0.15 |    1.124 |    1.1484 |       1.1490 |         1.0214 | 0.0680 |    0.0071 |  0.0908 |     0.0095 |
| ehr   | eur        |  310522 |     876421 |     0.26160 |     0.15 |    1.444 |    1.7548 |       1.5807 |         1.0279 | 0.0313 |    0.0011 |  0.0485 |     0.0018 |
| quest | eur        |   55519 |     346777 |     0.13800 |     0.15 |    1.194 |    1.2725 |       1.2266 |         1.0156 | 0.0328 |    0.0018 |  0.0824 |     0.0045 |
| self  | eur        |  114992 |    1789651 |     0.06037 |     0.15 |    1.453 |    1.9058 |       1.6485 |         0.9803 | 0.0240 |    0.0009 |  0.1265 |     0.0049 |

Genetic correlations

``` r
corrplot(cov2cor(mdd_covstruct$S), method='number')
```

![](/gpfs/igmmfs01/eddie/GenScotDepression/madams23/projects/mdd-meta/docs/gsem_files/figure-gfm/rg-1.png)<!-- -->

![](gsem_files/figure-gfm/rg-1.png)
