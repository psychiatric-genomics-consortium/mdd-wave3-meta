Conditional and Joint Analysis
================

``` r
library(rtracklayer)
library(readxl)
library(readr)
library(dplyr)
library(stringr)
```

# Methods

We ran a [conditional and joint
analysis](https://www.nature.com/articles/ng.2213) using
[GCTA](https://cnsgenomics.com/software/gcta/#COJO) to refine the list
of independent loci.

-   Final meta-analysed SNPs from the Ricopili pipeline were used.
-   Ricopili was used for initial clumping with index SNPs identified
    with *p*\< 10^{-4} and *r*<sup>2</sup>\< 0.1 within 3000kb windows.
    The extended MHC region was clumped as a single region.
-   Sumstats were filtered for MAF >= 0.01, INFO > 0.6, and Neff >= 80%
    of max.
-   Regions with a genome-wide significant SNP (p \< 5-e8) were
    identified from the clumped results. Regions within 50kb of each
    other were merged.
-   SNPs from these regions were extracted and filtered to unrelated of
    self- and genotype-identified European ancestry participants from UK
    Biobank.
-   A conditional analysis was performed on each region using the
    filtered sumstats superimposed on the UK Biobank LD structure.
-   Singleton regions (with only one genome-wide significant variant)
    were removed.

## Previous sumstats

Variants from [Howard et al
2019](https://www.nature.com/articles/s41593-018-%200326-7) (remove
first and last rows with table captions) and [Levey et al
2021](https://doi.org/10.1038/s41593-021-00860-2)

``` r
howard <- read_excel(snakemake@input$howard, skip=2, n_max=102)
```

    ## New names:
    ## * `Odds Ratio` -> `Odds Ratio...8`
    ## * `Lower 95% Confidence Interval` -> `Lower 95% Confidence Interval...9`
    ## * `Upper 95% Confidence Interval` -> `Upper 95% Confidence Interval...10`
    ## * `Log(Odds Ratio)` -> `Log(Odds Ratio)...11`
    ## * `Standard error of the Log(Odds Ratio)` -> `Standard error of the Log(Odds Ratio)...12`
    ## * ...

``` r
levey <- read_tsv(snakemake@input$levey, col_types=cols(CHR.BP=col_character()))
```

# Results

## COJO SNP and region counts

List of clumped and COJO SNPs and regions

-   Sumstats:
    results/cojo/daner_pgc_mdd_full_eur_hg19_v3.49.24.09.qc.gz  
-   Clump file:
    results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.09.gz.p4.clump.areator.sorted.1mhc  
-   COJO regions: 577  
-   Clumped SNPs: 817  
-   COJO Selected SNPs: 617  
-   Singleton regions: 65  
-   COJO+Clump SNPs: 617  
-   COJO Final SNPs: 552  
-   COJO Final SNPs p \<= 5e-8, pJ \<= 5e-8: 543  
-   COJO Final SNPs p > 5e-8, pJ \<= 5e-8: 9

Load list of COJO SNPs

``` r
cojo <- read_tsv(snakemake@input$cojo)
```

    ## 
    ## ── Column specification ───────────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   SNP = col_character(),
    ##   A1 = col_character(),
    ##   A2 = col_character(),
    ##   Direction = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

## Clumped results

Clumped results from Ricopili

``` r
rp <- read_table2(snakemake@input$rp_clump) %>% filter(P <= 5e-8)
```

    ## 
    ## ── Column specification ───────────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   SNP = col_character(),
    ##   A1A2 = col_character(),
    ##   `(Nca,Nco,Neff)Dir` = col_character(),
    ##   `LD-friends(0.1).p0.001` = col_character(),
    ##   `LD-friends(0.6).p0.001` = col_character(),
    ##   gwas_catalog_span.6 = col_character(),
    ##   `genes.6.50kb(dist2index)` = col_character(),
    ##   N.genes.6.50kb = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: 94 parsing failures.
    ##  row col expected actual                                                                                           file
    ## 1420 ngt a double      - 'results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.09.gz.p4.clump.areator.sorted.1mhc'
    ## 1850 ngt a double      - 'results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.09.gz.p4.clump.areator.sorted.1mhc'
    ## 1919 ngt a double      - 'results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.09.gz.p4.clump.areator.sorted.1mhc'
    ## 2244 ngt a double      - 'results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.09.gz.p4.clump.areator.sorted.1mhc'
    ## 2345 ngt a double      - 'results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.09.gz.p4.clump.areator.sorted.1mhc'
    ## .... ... ........ ...... ..............................................................................................
    ## See problems(...) for more details.

## Genomic ranges

Construct genomic range objects so that SNPs and regions can be
intersected and compared.

``` r
cojo_gr <- with(cojo, GRanges(seqnames=CHR, ranges=IRanges(start=range.left, end=range.right), SNP=SNP))
rp_gr <- with(rp, GRanges(seqnames=CHR, ranges=IRanges(start=range.left, end=range.right), SNP=SNP))
howard_gr <- with(howard, GRanges(seqnames=Chromosome, ranges=IRanges(start=`Postion (bp)`, width=1)))
levey_gr <- with(levey, GRanges(seqnames=CHR, ranges=IRanges(start=BP, width=1)))
```

## GWAS catalog

Tally entries in the GWAS catalog for each region

``` r
parse_catalog_entry <- function(catalog_entry) {
    # extract LD and RSID from first element
    ld_snp_match <- str_match(catalog_entry[1], "\\((.+)\\)(rs[[:digit:]]+)")
    ld <- as.numeric(ld_snp_match[,2])
    snp <- ld_snp_match[,3]
    # parse rest of elements into phenotype(PubMed ID)(P-value)
    phenotype_matches <- str_match(catalog_entry[-1], "(.+)\\(([[:digit:]]+)\\)\\((.+)\\)")
    phenotypes <- phenotype_matches[,2]
    pubmeds <- as.numeric(phenotype_matches[,3])
    p_values <- as.numeric(phenotype_matches[,4])
    
    return(data.frame(ld, catalogSNP=snp, phenotype=phenotypes, pubmed_id=pubmeds, P=p_values))
}


rp_gwas_catalog_entries <-
plyr::ddply(filter(rp, gwas_catalog_span.6 != '-'), ~SNP, function(rp_entry) {
    # get GWAS catalog cell
    gwas_catalog <- rp_entry$gwas_catalog_span.6
    # split out into catalog entries (separated by /, removing first empty element)
    catalog_entries <- str_split(gwas_catalog, "/")[[1]]
    catalog_entries_complete <- catalog_entries[which(catalog_entries != "")]
    # split SNP entries by ";"
    catalog_entries_list <- str_split(catalog_entries_complete, ";")
    # parse each entry
    catalog_entries_df <- plyr::ldply(catalog_entries_list, parse_catalog_entry)
    return(catalog_entries_df)
}) %>% as_tibble()
```

# Parse gene list

``` r
rp_genes_dist <- 
plyr::ddply(filter(rp, `genes.6.50kb(dist2index)` != '-'), ~SNP, function(rp_entry) {
    genes_dist <- rp_entry$`genes.6.50kb(dist2index)`
    genes_dist_list <- str_split(genes_dist, ',')[[1]]
    genes_dist_match <- str_match(genes_dist_list, "(.+)\\((.+)\\)")
    return(data.frame(gene=genes_dist_match[,2], dist2index=as.numeric(genes_dist_match[,3])))
}) %>% as_tibble()
```

## Comparison to previous findings

Find which COJO regions overlap with Howard

``` r
cojo_howard_overlaps <- findOverlaps(cojo_gr, howard_gr)
cojo_howard_overlaps
```

    ## Hits object with 115 hits and 0 metadata columns:
    ##         queryHits subjectHits
    ##         <integer>   <integer>
    ##     [1]         1           1
    ##     [2]         5           2
    ##     [3]         6           3
    ##     [4]         8           4
    ##     [5]         9           4
    ##     ...       ...         ...
    ##   [111]       521          97
    ##   [112]       525          98
    ##   [113]       532         100
    ##   [114]       534         101
    ##   [115]       548         102
    ##   -------
    ##   queryLength: 552 / subjectLength: 102

Count number of regions in Howard that overlap:

``` r
howard %>% slice(unique(cojo_howard_overlaps@to)) %>% count()
```

    ## # A tibble: 1 × 1
    ##       n
    ##   <int>
    ## 1    88

Find which COJO regions overlap with Levey

``` r
cojo_levey_overlaps <- findOverlaps(cojo_gr, levey_gr)
cojo_levey_overlaps
```

    ## Hits object with 260 hits and 0 metadata columns:
    ##         queryHits subjectHits
    ##         <integer>   <integer>
    ##     [1]         1         110
    ##     [2]         5           2
    ##     [3]         6          65
    ##     [4]         8          77
    ##     [5]         8          26
    ##     ...       ...         ...
    ##   [256]       537         154
    ##   [257]       538         154
    ##   [258]       544         120
    ##   [259]       548         153
    ##   [260]       549         211
    ##   -------
    ##   queryLength: 552 / subjectLength: 223

Count number of regions in Levey that overlap:

``` r
levey %>% slice(unique(cojo_levey_overlaps@to)) %>% count()
```

    ## # A tibble: 1 × 1
    ##       n
    ##   <int>
    ## 1   191

``` r
cojo_known <- cojo %>% slice(unique(cojo_levey_overlaps@from))

catalog_known <- rp_gwas_catalog_entries %>% filter(SNP %in% cojo_known$SNP) %>% count(phenotype) %>% arrange(desc(n))
catalog_known
```

    ## # A tibble: 166 × 2
    ##    phenotype                   n
    ##    <chr>                   <int>
    ##  1 Serum_metabolit...        186
    ##  2 Schizophrenia              45
    ##  3 Intelligence_(MTAG)        42
    ##  4 Trans_fatty_acid_levels    29
    ##  5 Waist_circumfer...         29
    ##  6 Depression_(broad)         28
    ##  7 Glycerophosphol...         27
    ##  8 Neuroticism                26
    ##  9 Autism_spectrum...         24
    ## 10 Depression                 22
    ## # … with 156 more rows

Newly discovered regions

``` r
cojo_new <- cojo %>% slice(unique(-cojo_levey_overlaps@from)) %>% arrange(P)
cojo_new %>%
select(region, snp_idx, CHR, SNP, BP, P, pJ) %>%
group_by(region)
```

    ## # A tibble: 360 × 7
    ## # Groups:   region [342]
    ##    region snp_idx   CHR SNP               BP        P       pJ
    ##     <dbl>   <dbl> <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ##  1    174       1     5 rs1993739  153215007 5.71e-20 5.72e-20
    ##  2    487       1    20 rs17805843  51205902 9.34e-20 9.36e-20
    ##  3     49       2     2 rs359247    60477052 2.67e-18 6.60e-17
    ##  4    166       1     5 rs6863440  120084621 3.25e-18 3.26e-18
    ##  5    211       1     7 rs10499337   3521803 1.19e-17 1.48e-16
    ##  6    362       1    12 rs2363585   60791165 1.69e-16 1.69e-16
    ##  7    310       1    10 rs12778915  77617557 2.35e-16 2.35e-16
    ##  8    374       1    12 rs7962128  121907336 3.68e-16 3.68e-16
    ##  9     76       1     2 rs13418032 198413692 9.89e-16 9.90e-16
    ## 10     81       1     2 rs72931605 212693775 1.29e-15 1.21e- 9
    ## # … with 350 more rows

``` r
rp_gwas_catalog_entries %>% filter(SNP %in% cojo_new$SNP) %>% count(phenotype) %>% arrange(desc(n)) %>% filter(!phenotype %in% catalog_known$phenotype)
```

    ## # A tibble: 61 × 2
    ##    phenotype                    n
    ##    <chr>                    <int>
    ##  1 Body_mass_index...           9
    ##  2 Morning_vs._eve...           5
    ##  3 Obesity                      5
    ##  4 Hip_circumference            4
    ##  5 Diastolic_blood_pressure     3
    ##  6 Monocyte_count               3
    ##  7 Pulse_pressure               3
    ##  8 Ulcerative_colitis           3
    ##  9 Waist_circumference          3
    ## 10 Adiponectin_levels           2
    ## # … with 51 more rows

``` r
rp_genes_dist %>% filter(SNP %in% cojo_new$SNP) %>% group_by(SNP) %>% filter(dist2index == min(dist2index)) %>% ungroup() %>% select(gene) %>% distinct()
```

    ## # A tibble: 374 × 1
    ##    gene    
    ##    <chr>   
    ##  1 NISCH   
    ##  2 STAB1   
    ##  3 NT5DC2  
    ##  4 SMIM4   
    ##  5 PBRM1   
    ##  6 NTRK3   
    ##  7 GTF2IRD1
    ##  8 NXPH1   
    ##  9 MAGI2   
    ## 10 ETV6    
    ## # … with 364 more rows

## Clumped results

Find regions overlapping between clumped and COJO results

``` r
cojo_clumped_overlaps <- findOverlaps(cojo_gr, rp_gr)
```

List COJO results where the selected SNP is not in the clumped results

``` r
cojo_newly_selected <- 
cojo %>% slice(unique(cojo_clumped_overlaps@from)) %>%
filter(!SNP %in% rp$SNP) %>%
select(region, snp_idx, CHR, SNP, P, pJ) 
cojo_newly_selected
```

    ## # A tibble: 58 × 6
    ##    region snp_idx   CHR SNP                P       pJ
    ##     <dbl>   <dbl> <dbl> <chr>          <dbl>    <dbl>
    ##  1      1       1     1 rs301806    1.87e-16 1.87e-16
    ##  2     10       1     1 rs437021    5.70e-11 5.71e-11
    ##  3     13       1     1 rs3101341   5.22e-27 2.29e-11
    ##  4     13       2     1 rs2797104   5.84e-32 4.16e-10
    ##  5     47       1     2 rs858938    3.87e- 8 2.65e-10
    ##  6     48       2     2 rs56873970  2.76e- 3 1.82e- 8
    ##  7     49       1     2 rs7576017   1.18e-13 2.96e-12
    ##  8     50       1     2 rs149044563 7.75e-13 7.76e-13
    ##  9     53       1     2 rs73949838  5.43e- 9 5.43e- 9
    ## 10     55       1     2 rs3860446   1.10e-10 3.29e-11
    ## # … with 48 more rows

Newly selected SNPs that were not GWsig in the clumped results

``` r
cojo_newly_selected %>%
filter(P > 5e-8)
```

    ## # A tibble: 9 × 6
    ##   region snp_idx   CHR SNP                   P       pJ
    ##    <dbl>   <dbl> <dbl> <chr>             <dbl>    <dbl>
    ## 1     48       2     2 rs56873970 0.00276      1.82e- 8
    ## 2     60       1     2 rs2381462  0.00000896   9.62e-11
    ## 3     69       2     2 rs1371187  0.0000000835 4.50e-11
    ## 4    118       1     3 rs13073224 0.0000000573 2.02e- 8
    ## 5    186       2     6 rs2747467  0.144        1.63e-12
    ## 6    248       1     8 rs10503484 0.000000367  7.61e- 9
    ## 7    341       4    11 rs73004019 0.000154     5.15e-10
    ## 8    386       2    13 rs2329076  0.000000154  3.24e- 8
    ## 9    412       2    15 rs4774501  0.00000214   6.16e-11

List clumped SNPs in retained regions that were not selected by COJO

``` r
rp %>% slice(unique(cojo_clumped_overlaps@to)) %>%
filter(!SNP %in% cojo$SNP)
```

    ## # A tibble: 246 × 23
    ##    SNP          CHR       BP        P    OR     SE A1A2  FRQ_A_524857 FRQ_U_3059006
    ##    <chr>      <dbl>    <dbl>    <dbl> <dbl>  <dbl> <chr>        <dbl>         <dbl>
    ##  1 rs301817       1  8503379 3.23e-17 1.02  0.0028 C/A         0.419         0.426 
    ##  2 rs75986133     1 52846058 1.99e- 8 0.967 0.0061 A/G         0.044         0.0458
    ##  3 rs446952       1 61738636 1.78e-11 1.02  0.0027 T/C         0.465         0.458 
    ##  4 rs2568957      1 72764430 1.58e-32 0.963 0.0031 A/G         0.176         0.172 
    ##  5 rs12127789     1 72740073 4.91e-23 1.04  0.0038 T/G         0.119         0.119 
    ##  6 rs12748090     1 73005277 1.47e-14 1.02  0.003  A/T         0.279         0.277 
    ##  7 rs75805282     1 73857518 6.32e-13 0.950 0.0071 T/C         0.0301        0.0304
    ##  8 rs61771936     1 73266056 3.66e-12 0.977 0.0034 A/G         0.156         0.158 
    ##  9 rs12128239     1 73059446 4.20e-12 0.977 0.0034 G/A         0.163         0.164 
    ## 10 rs10736420     1 74000649 1.18e-11 1.03  0.0036 C/T         0.12          0.118 
    ## # … with 236 more rows, and 14 more variables: INFO <dbl>,
    ## #   (Nca,Nco,Neff)Dir <chr>, ngt <dbl>, LD-friends(0.1).p0.001 <chr>,
    ## #   range.left <dbl>, range.right <dbl>, span(kb) <dbl>,
    ## #   LD-friends(0.6).p0.001 <chr>, range.left.6 <dbl>, range.right.6 <dbl>,
    ## #   span.6(kb) <dbl>, gwas_catalog_span.6 <chr>,
    ## #   genes.6.50kb(dist2index) <chr>, N.genes.6.50kb <chr>

Line up non-selected SNPs with selected SNPs in the region

``` r
bind_cols(
select(slice(cojo, cojo_clumped_overlaps@from), region, snp_idx, SNP.cojo=SNP, BP.cojo=BP, P.cojo=P, PJ.cojo=pJ),
select(slice(rp, cojo_clumped_overlaps@to), SNP.rp=SNP, BP.rp=BP, P.rp=P)
) %>%
filter(!SNP.rp %in% cojo$SNP)
```

    ## # A tibble: 400 × 9
    ##    region snp_idx SNP.cojo   BP.cojo   P.cojo  PJ.cojo SNP.rp      BP.rp     P.rp
    ##     <dbl>   <dbl> <chr>        <dbl>    <dbl>    <dbl> <chr>       <dbl>    <dbl>
    ##  1      1       1 rs301806   8482078 1.87e-16 1.87e-16 rs301817   8.50e6 3.23e-17
    ##  2      9       1 rs7413471 52339759 2.96e-15 2.96e-15 rs75986133 5.28e7 1.99e- 8
    ##  3     10       1 rs437021  61738270 5.70e-11 5.71e-11 rs446952   6.17e7 1.78e-11
    ##  4     13       1 rs3101341 72747844 5.22e-27 2.29e-11 rs2568957  7.28e7 1.58e-32
    ##  5     13       1 rs3101341 72747844 5.22e-27 2.29e-11 rs12127789 7.27e7 4.91e-23
    ##  6     13       1 rs3101341 72747844 5.22e-27 2.29e-11 rs12748090 7.30e7 1.47e-14
    ##  7     13       1 rs3101341 72747844 5.22e-27 2.29e-11 rs75805282 7.39e7 6.32e-13
    ##  8     13       1 rs3101341 72747844 5.22e-27 2.29e-11 rs61771936 7.33e7 3.66e-12
    ##  9     13       1 rs3101341 72747844 5.22e-27 2.29e-11 rs12128239 7.31e7 4.20e-12
    ## 10     13       1 rs3101341 72747844 5.22e-27 2.29e-11 rs10736420 7.40e7 1.18e-11
    ## # … with 390 more rows
