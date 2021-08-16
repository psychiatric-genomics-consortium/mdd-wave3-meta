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
-   Sumstats were filtered for MAF >= 0.01 and INFO > 0.6
-   Regions with a genome-wide significant SNP (p \< 5-e8) were
    identified from the clumped results. Regions within 50kb of each
    other were merged.
-   SNPs from these regions were extracted filtered to unrelated of
    self- and genotype-identified European ancestry participants from UK
    Biobank.
-   A conditional analysis was performed on each region using the
    filtered sumstats superimposed on the UK Biobank LD structure.

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
    results/cojo/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp.qc.gz  
-   Clump file:
    results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp.gz.p4.clump.areator.sorted.1mhc  
-   COJO regions: 549  
-   Clumped SNPs: 738  
-   COJO Selected SNPs: 594  
-   Singleton regions: 63  
-   COJO+Clump SNPs: 606  
-   COJO Final SNPs: 543  
-   COJO Final SNPs p \<= 5e-8, pJ \<= 5e-8: 512  
-   COJO Final SNPs p > 5e-8, pJ \<= 5e-8: 29

Load list of COJO SNPs

``` r
cojo <- read_tsv(snakemake@input$cojo)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   SNP = col_character(),
    ##   A1 = col_character(),
    ##   A2 = col_character(),
    ##   Direction = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

## Clumped results

Clumped results from Ricopili and METACARPA for comparison

``` r
rp <- read_table2(snakemake@input$rp_clump) %>% filter(P <= 5e-8)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   SNP = col_character(),
    ##   A1A2 = col_character(),
    ##   `(Nca,Nco,Neff)Dir` = col_character(),
    ##   ngt = col_character(),
    ##   `LD-friends(0.1).p0.001` = col_character(),
    ##   `LD-friends(0.6).p0.001` = col_character(),
    ##   gwas_catalog_span.6 = col_character(),
    ##   `genes.6.50kb(dist2index)` = col_character(),
    ##   N.genes.6.50kb = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
mc <- read_table2(snakemake@input$mc_clump) %>% filter(P <= 5e-8)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   SNP = col_character(),
    ##   A1A2 = col_character(),
    ##   ngt = col_character(),
    ##   `LD-friends(0.1).p0.001` = col_character(),
    ##   `LD-friends(0.6).p0.001` = col_character(),
    ##   gwas_catalog_span.6 = col_character(),
    ##   `genes.6.50kb(dist2index)` = col_character(),
    ##   N.genes.6.50kb = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

## Genomic ranges

Construct genomic range objects so that SNPs and regions can be
intersected and compared.

``` r
cojo_gr <- with(cojo, GRanges(seqnames=CHR, ranges=IRanges(start=range.left, end=range.right), SNP=SNP))
rp_gr <- with(rp, GRanges(seqnames=CHR, ranges=IRanges(start=range.left, end=range.right), SNP=SNP))
mc_gr <- with(mc, GRanges(seqnames=CHR, ranges=IRanges(start=range.left, end=range.right), SNP=SNP))
howard_gr <- with(howard, GRanges(seqnames=Chromosome, ranges=IRanges(start=`Postion (bp)`, width=1)))
levey_gr <- with(levey, GRanges(seqnames=CHR, ranges=IRanges(start=BP, width=1)))
```

``` r
rp_mc_overlaps <- findOverlaps(rp_gr, mc_gr)

rp_and_mc <- rp %>% slice(unique(rp_mc_overlaps@from))
rp_not_mc <- rp %>% slice(-unique(rp_mc_overlaps@from))
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

    ## Hits object with 128 hits and 0 metadata columns:
    ##         queryHits subjectHits
    ##         <integer>   <integer>
    ##     [1]         1           1
    ##     [2]         5           2
    ##     [3]         6           3
    ##     [4]         9           4
    ##     [5]        10           4
    ##     ...       ...         ...
    ##   [124]       519          97
    ##   [125]       522          98
    ##   [126]       527         100
    ##   [127]       528         101
    ##   [128]       541         102
    ##   -------
    ##   queryLength: 543 / subjectLength: 102

Count number of regions in Howard that overlap:

``` r
howard %>% slice(unique(cojo_howard_overlaps@to)) %>% count()
```

    ## # A tibble: 1 × 1
    ##       n
    ##   <int>
    ## 1    85

Find which COJO regions overlap with Levey

``` r
cojo_levey_overlaps <- findOverlaps(cojo_gr, levey_gr)
cojo_levey_overlaps
```

    ## Hits object with 270 hits and 0 metadata columns:
    ##         queryHits subjectHits
    ##         <integer>   <integer>
    ##     [1]         1         110
    ##     [2]         5           2
    ##     [3]         6          65
    ##     [4]         9          77
    ##     [5]         9          26
    ##     ...       ...         ...
    ##   [266]       531         154
    ##   [267]       532         154
    ##   [268]       538         120
    ##   [269]       541         153
    ##   [270]       542         211
    ##   -------
    ##   queryLength: 543 / subjectLength: 223

Count number of regions in Levey that overlap:

``` r
levey %>% slice(unique(cojo_levey_overlaps@to)) %>% count()
```

    ## # A tibble: 1 × 1
    ##       n
    ##   <int>
    ## 1   186

``` r
cojo_known <- cojo %>% slice(unique(cojo_levey_overlaps@from))

catalog_known <- rp_gwas_catalog_entries %>% filter(SNP %in% cojo_known$SNP) %>% count(phenotype) %>% arrange(desc(n))
catalog_known
```

    ## # A tibble: 76 × 2
    ##    phenotype                     n
    ##    <chr>                     <int>
    ##  1 Schizophrenia                26
    ##  2 Neuroticism                  18
    ##  3 Depression                   17
    ##  4 Intelligence_(MTAG)          16
    ##  5 Depressive_symp...           15
    ##  6 Body_mass_index              12
    ##  7 Autism_spectrum...           11
    ##  8 Depressive_symptoms          10
    ##  9 Neuroticism_(MTAG)           10
    ## 10 Major_depressive_disorder     9
    ## # … with 66 more rows

Newly discovered regions

``` r
cojo_new <- cojo %>% slice(unique(-cojo_levey_overlaps@from)) %>% arrange(P)
cojo_new
```

    ## # A tibble: 351 × 28
    ##    region snp_idx   CHR SNP               BP A1    A2    FRQ_A_511755 FRQ_U_3072803
    ##     <dbl>   <dbl> <dbl> <chr>          <dbl> <chr> <chr>        <dbl>         <dbl>
    ##  1    476       1    20 rs34966255  51193862 C     T            0.18          0.183
    ##  2     50       2     2 rs951807    60513389 T     C            0.404         0.405
    ##  3    182       1     5 rs1993739  153215007 T     C            0.166         0.171
    ##  4     50       1     2 rs10490064  60163127 A     G            0.367         0.38 
    ##  5    215       1     7 rs1114780    3551788 C     T            0.237         0.242
    ##  6    243       1     7 rs7797141  115052164 A     G            0.274         0.266
    ##  7    174       1     5 rs7726836  120065919 T     A            0.294         0.285
    ##  8    354       1    12 rs2363585   60791165 A     T            0.375         0.383
    ##  9    368       1    12 rs7962128  121907336 A     G            0.455         0.46 
    ## 10     80       1     2 rs13418032 198413692 G     C            0.334         0.323
    ## # … with 341 more rows, and 19 more variables: INFO <dbl>, OR <dbl>, SE <dbl>,
    ## #   P <dbl>, ngt <dbl>, Direction <chr>, HetISqt <dbl>, HetDf <dbl>,
    ## #   HetPVa <dbl>, Nca <dbl>, Nco <dbl>, Neff_half <dbl>, bJ <dbl>, bJ_se <dbl>,
    ## #   pJ <dbl>, LD_r <dbl>, range.left <dbl>, range.right <dbl>, n_sig_snps <dbl>

``` r
rp_gwas_catalog_entries %>% filter(SNP %in% cojo_new$SNP) %>% count(phenotype) %>% arrange(desc(n)) %>% filter(!phenotype %in% catalog_known$phenotype)
```

    ## # A tibble: 52 × 2
    ##    phenotype                    n
    ##    <chr>                    <int>
    ##  1 Height                       4
    ##  2 Diastolic_blood_pressure     3
    ##  3 Morning_vs._eve...           3
    ##  4 Adiponectin_levels           2
    ##  5 Blood_osmolalit...           2
    ##  6 Hip_circumferen...           2
    ##  7 LDL_cholesterol_levels       2
    ##  8 Lobe_attachment...           2
    ##  9 Lymphocyte_perc...           2
    ## 10 Mean_corpuscular_volume      2
    ## # … with 42 more rows

``` r
rp_genes_dist %>% filter(SNP %in% cojo_new$SNP) %>% group_by(SNP) %>% filter(dist2index == min(dist2index)) %>% ungroup() %>% select(gene) %>% distinct()
```

    ## # A tibble: 261 × 1
    ##    gene    
    ##    <chr>   
    ##  1 PPP3CA  
    ##  2 MIR1255A
    ##  3 FLJ20021
    ##  4 NISCH   
    ##  5 STAB1   
    ##  6 NT5DC2  
    ##  7 SMIM4   
    ##  8 PBRM1   
    ##  9 NTRK3   
    ## 10 GTF2IRD1
    ## # … with 251 more rows

## Comparison between pre-COJO and post-COJO

``` r
rp_genes_dist %>% filter(!SNP %in% cojo$SNP)
```

    ## # A tibble: 1,108 × 3
    ##    SNP        gene       dist2index
    ##    <chr>      <chr>           <dbl>
    ##  1 rs10064584 LINC02240           0
    ##  2 rs10089596 SGCZ                0
    ##  3 rs10165927 PLA2R1              0
    ##  4 rs10165927 ITGB6               0
    ##  5 rs10170048 UBE2F-SCLY          0
    ##  6 rs10170048 SCLY                0
    ##  7 rs10170048 ESPNL               0
    ##  8 rs10170048 KLHL30              0
    ##  9 rs10170048 ERFE                0
    ## 10 rs10170048 ILKAP               0
    ## # … with 1,098 more rows

## Comparison between Ricopili and METACARPA clumped results

Shared regions and regions not significant in the METACARPA analysis.

``` r
rp_and_mc_phenotypes <- rp_gwas_catalog_entries %>% filter(SNP %in% rp_and_mc$SNP) %>% count(phenotype) %>% arrange(desc(n))
rp_not_mc_phenotypes <- rp_gwas_catalog_entries %>% filter(SNP %in% rp_not_mc$SNP) %>% count(phenotype) %>% arrange(desc(n))

rp_not_mc_phenotypes %>% filter(!phenotype %in% rp_and_mc_phenotypes$phenotype)
```

    ## # A tibble: 10 × 2
    ##    phenotype                    n
    ##    <chr>                    <int>
    ##  1 Menopause_(age_at_onset)     3
    ##  2 Mean_corpuscular_volume      2
    ##  3 Visceral_adipos...           2
    ##  4 Adenocarcinoma               1
    ##  5 Hemoglobin                   1
    ##  6 Lipoprotein-ass...           1
    ##  7 Mean_corpuscula...           1
    ##  8 Monobrow                     1
    ##  9 Prostate_cancer              1
    ## 10 White_blood_cel...           1

``` r
rp_genes_dist %>% group_by(SNP) %>% filter(dist2index == min(dist2index)) %>% ungroup() %>% select(gene) %>% distinct()
```

    ## # A tibble: 790 × 1
    ##    gene     
    ##    <chr>    
    ##  1 PPP3CA   
    ##  2 MIR1255A 
    ##  3 FLJ20021 
    ##  4 CD14     
    ##  5 LINC02240
    ##  6 SGCZ     
    ##  7 NISCH    
    ##  8 STAB1    
    ##  9 NT5DC2   
    ## 10 SMIM4    
    ## # … with 780 more rows

``` r
rp_genes_dist %>% filter(SNP %in% rp_and_mc$SNP) %>% group_by(SNP) %>% filter(dist2index == min(dist2index)) %>% ungroup() %>% select(gene) %>% distinct()
```

    ## # A tibble: 679 × 1
    ##    gene     
    ##    <chr>    
    ##  1 PPP3CA   
    ##  2 MIR1255A 
    ##  3 FLJ20021 
    ##  4 CD14     
    ##  5 LINC02240
    ##  6 SGCZ     
    ##  7 NISCH    
    ##  8 STAB1    
    ##  9 NT5DC2   
    ## 10 SMIM4    
    ## # … with 669 more rows

``` r
rp_genes_dist %>% filter(SNP %in% rp_not_mc$SNP) %>% group_by(SNP) %>% filter(dist2index == min(dist2index)) %>% ungroup() %>% select(gene) %>% distinct()
```

    ## # A tibble: 118 × 1
    ##    gene      
    ##    <chr>     
    ##  1 UBE2F-SCLY
    ##  2 SCLY      
    ##  3 ESPNL     
    ##  4 KLHL30    
    ##  5 ERFE      
    ##  6 ILKAP     
    ##  7 SHTN1     
    ##  8 VAX1      
    ##  9 MIR3663HG 
    ## 10 MIR3663   
    ## # … with 108 more rows
