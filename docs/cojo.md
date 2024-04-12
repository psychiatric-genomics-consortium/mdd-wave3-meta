Conditional and Joint Analysis
================

``` r
library(plyranges)
library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(UpSetR)
library(genpwr)
library(ggplot2)
```

# Methods

We ran a [conditional and joint
analysis](https://www.nature.com/articles/ng.2213) using
[GCTA](https://cnsgenomics.com/software/gcta/#COJO) to refine the list
of independent loci.

- Final meta-analysed SNPs from the Ricopili pipeline were used.
- Ricopili was used for initial clumping with index SNPs identified with
  $p <$ 10^{-4} and $r^2 <$ 0.1 within 3000kb windows. The extended MHC
  region was clumped as a single region.
- Sumstats were filtered for MAF \>= 0.01, INFO \> 0.6, and Neff \>= 80%
  of max.
- Regions with a genome-wide significant SNP (p \< 5-e8) were identified
  from the clumped results.
- SNPs from these regions were extracted and filtered to 7500 unrelated
  of self- and genotype-identified European ancestry participants from
  UK Biobank.
- A conditional analysis was performed on with 10MB windows using the
  filtered sumstats superimposed on the UK Biobank LD structure.

## Previous sumstats

Variants from [Wray et al
2018](https://www.nature.com/articles/s41588-018-0090-3/tables/2),
[Howard et al
2019](https://www.nature.com/articles/s41593-018-%200326-7) (remove
first and last rows with table captions), [Levey et al
2021](https://doi.org/10.1038/s41593-021-00860-2), [Giannakopoulou et al
2021](https://jamanetwork.com/journals/jamapsychiatry/article-abstract/2784695),
[Als et al
2022](https://www.medrxiv.org/content/10.1101/2022.08.24.22279149v1),
[Meng et al 2024](https://www.nature.com/articles/s41588-023-01596-4),
and the [GWAS catalog for unipolar
depression](https://www.ebi.ac.uk/gwas/efotraits/EFO_0003761). Parse out
regions from Wray and Howard and load queried regions for other results

``` r
tags <- read_table(snakemake@input$tags)
```

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   SNP = col_character(),
    ##   CHR = col_double(),
    ##   BP = col_double(),
    ##   NTAG = col_double(),
    ##   LEFT = col_double(),
    ##   RIGHT = col_double(),
    ##   KBSPAN = col_double(),
    ##   TAGS = col_character()
    ## )

``` r
wray <- read_tsv(snakemake@input$wray) %>%
    separate(`Region (Mb)`, into=c('range.left.Mb', 'range.right.Mb'), sep='–', convert=TRUE) %>%
    mutate(range.left=range.left.Mb*1e6, range.right=range.right.Mb*1e6)
```

    ## Rows: 44 Columns: 11

    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): Region (Mb), SNP, P, A1/A2, Prev., Gene context
    ## dbl (4): Chr., OR (A1), s.e. (log(OR)), Freq.
    ## num (1): Location (bp)
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
howard <- read_excel(snakemake@input$howard, skip=2, n_max=102) %>%
    separate(`Region (bp) of clumped variants (P < 10-4) in linkage disequilibrium (r2 > 0.1) with lead variant`, into=c('range.left', 'range.right'), sep='-', convert=TRUE)
```

    ## New names:
    ## • `Odds Ratio` -> `Odds Ratio...8`
    ## • `Lower 95% Confidence Interval` -> `Lower 95% Confidence Interval...9`
    ## • `Upper 95% Confidence Interval` -> `Upper 95% Confidence Interval...10`
    ## • `Log(Odds Ratio)` -> `Log(Odds Ratio)...11`
    ## • `Standard error of the Log(Odds Ratio)` -> `Standard error of the Log(Odds Ratio)...12`
    ## • `P-value` -> `P-value...13`
    ## • `Odds Ratio` -> `Odds Ratio...14`
    ## • `Lower 95% Confidence Interval` -> `Lower 95% Confidence Interval...15`
    ## • `Upper 95% Confidence Interval` -> `Upper 95% Confidence Interval...16`
    ## • `Log(Odds Ratio)` -> `Log(Odds Ratio)...17`
    ## • `Standard error of the Log(Odds Ratio)` -> `Standard error of the Log(Odds Ratio)...18`
    ## • `P-value` -> `P-value...19`
    ## • `Odds Ratio` -> `Odds Ratio...20`
    ## • `Lower 95% Confidence Interval` -> `Lower 95% Confidence Interval...21`
    ## • `Upper 95% Confidence Interval` -> `Upper 95% Confidence Interval...22`
    ## • `Log(Odds Ratio)` -> `Log(Odds Ratio)...23`
    ## • `Standard error of the Log(Odds Ratio)` -> `Standard error of the Log(Odds Ratio)...24`
    ## • `P-value` -> `P-value...25`
    ## • `Odds Ratio` -> `Odds Ratio...26`
    ## • `Lower 95% Confidence Interval` -> `Lower 95% Confidence Interval...27`
    ## • `Upper 95% Confidence Interval` -> `Upper 95% Confidence Interval...28`
    ## • `Log(Odds Ratio)` -> `Log(Odds Ratio)...29`
    ## • `Standard error of the Log(Odds Ratio)` -> `Standard error of the Log(Odds Ratio)...30`
    ## • `P-value` -> `P-value...31`
    ## • `Direction` -> `Direction...32`
    ## • `Odds Ratio` -> `Odds Ratio...33`
    ## • `Lower 95% Confidence Interval` -> `Lower 95% Confidence Interval...34`
    ## • `Upper 95% Confidence Interval` -> `Upper 95% Confidence Interval...35`
    ## • `Log(Odds Ratio)` -> `Log(Odds Ratio)...36`
    ## • `Standard error of the Log(Odds Ratio)` -> `Standard error of the Log(Odds Ratio)...37`
    ## • `P-value` -> `P-value...38`
    ## • `Odds Ratio` -> `Odds Ratio...39`
    ## • `Lower 95% Confidence Interval` -> `Lower 95% Confidence Interval...40`
    ## • `Upper 95% Confidence Interval` -> `Upper 95% Confidence Interval...41`
    ## • `Log(Odds Ratio)` -> `Log(Odds Ratio)...42`
    ## • `Standard error of the Log(Odds Ratio)` -> `Standard error of the Log(Odds Ratio)...43`
    ## • `P-value` -> `P-value...44`
    ## • `Direction` -> `Direction...45`

``` r
levey <- read_tsv(snakemake@input$levey, col_types=cols(CHR.BP=col_character())) %>%
left_join(select(tags, SNP, LEFT, RIGHT), by=c('rsid'='SNP')) %>%
mutate(LEFT=if_else(is.na(LEFT), true=BP, false=LEFT),
       RIGHT=if_else(is.na(RIGHT), true=BP, false=RIGHT))

giannakopoulou <- read_tsv(snakemake@input$giannakopoulou, col_types=cols('CHR:POS'=col_character())) %>%
separate(`CHR:POS`, into=c('CHR', 'POS'), convert=TRUE) %>%
left_join(select(tags, SNP, LEFT, RIGHT), by='SNP') %>%
mutate(LEFT=if_else(is.na(LEFT), true=as.numeric(POS), false=LEFT),
       RIGHT=if_else(is.na(RIGHT), true=as.numeric(POS), false=RIGHT))
       
als <- read_excel(snakemake@input$als, sheet=1)

meng <- read_excel(snakemake@input$meng, sheet = 7, skip = 2) |>
  filter(!is.na(Chr)) |>
  select(SNP, CHR = Chr, BP = bp) |>
  left_join(select(tags, SNP, LEFT, RIGHT), by='SNP') |>
  mutate(LEFT = coalesce(LEFT, BP), RIGHT = coalesce(RIGHT, BP))
```

    ## New names:
    ## • `` -> `...6`
    ## • `` -> `...7`
    ## • `` -> `...8`
    ## • `` -> `...9`
    ## • `` -> `...10`
    ## • `` -> `...12`
    ## • `` -> `...13`
    ## • `` -> `...14`
    ## • `` -> `...15`
    ## • `` -> `...16`
    ## • `` -> `...18`
    ## • `` -> `...19`
    ## • `` -> `...20`
    ## • `` -> `...21`
    ## • `` -> `...22`
    ## • `` -> `...24`
    ## • `` -> `...25`
    ## • `` -> `...26`
    ## • `` -> `...27`
    ## • `` -> `...28`
    ## • `` -> `...30`
    ## • `` -> `...31`
    ## • `` -> `...32`
    ## • `` -> `...33`
    ## • `` -> `...34`
    ## • `` -> `...36`
    ## • `` -> `...37`
    ## • `` -> `...38`
    ## • `` -> `...39`
    ## • `` -> `...40`

``` r
gwas_catalog <- read_tsv(snakemake@input$catalog) %>%
filter(!is.na(CHR_ID)) %>%
mutate(CHR=if_else(CHR_ID == 'X', true=23, false=as.numeric(CHR_ID)),
       POS=as.numeric(CHR_POS),
       SNP=paste0('rs', SNP_ID_CURRENT)) %>%
filter(!is.na(POS)) %>%
left_join(select(tags, SNP, LEFT, RIGHT), by='SNP') %>%
mutate(LEFT=if_else(is.na(LEFT), true=POS, false=LEFT),
       RIGHT=if_else(is.na(RIGHT), true=POS, false=RIGHT))
```

    ## Rows: 2391 Columns: 38
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (28): FIRST AUTHOR, JOURNAL, LINK, STUDY, DISEASE/TRAIT, INITIAL SAMPLE...
    ## dbl   (8): PUBMEDID, UPSTREAM_GENE_DISTANCE, DOWNSTREAM_GENE_DISTANCE, MERGE...
    ## date  (2): DATE ADDED TO CATALOG, DATE
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## Warning: There were 2 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `CHR = if_else(CHR_ID == "X", true = 23, false = as.numeric(CHR_ID))`.
    ## Caused by warning in `if_else()`:
    ## ! NAs introduced by coercion
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.

# Results

## COJO SNP and region counts

Load list of COJO SNPs (with multiple or single GW-sig SNPs)

``` r
cojo <- read_tsv(snakemake@input$cojo)
```

    ## Rows: 697 Columns: 25
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (3): SNP, A1, A2
    ## dbl (22): Locus, SNP.Index, CHR, BP, FRQ_A, FRQ_U, INFO, OR, SE, P, HetISqt,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

## Clumped results

Clumped results from Ricopili

``` r
rp <- read_table(snakemake@input$rp_clump) %>% filter(P <= 5e-8)
```

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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

## Genomic ranges

Construct genomic range objects so that SNPs and regions can be
intersected and compared.

``` r
cojo_gr <- with(cojo, GRanges(seqnames=CHR, ranges=IRanges(start=Locus.Left, end=Locus.Right), SNP=SNP))

wray_gr <- with(wray, GRanges(seqnames=Chr., ranges=IRanges(start=range.left, end=range.right, SNP=SNP)))
howard_gr <- with(howard, GRanges(seqnames=Chromosome, ranges=IRanges(start=range.left, end=range.right, SNP=`Marker Name`)))
levey_gr <- with(levey, GRanges(seqnames=CHR, ranges=IRanges(start=LEFT, end=RIGHT)))
giannakopoulou_gr <- with(giannakopoulou, GRanges(seqnames=CHR, ranges=IRanges(start=LEFT, end=RIGHT)))
als_gr <- with(als, GRanges(seqnames=CHR, ranges=IRanges(start=range.left.index.SNP, end=range.right.index.SNP)))
meng_gr <- with(meng, GRanges(seqnames=CHR, ranges=IRanges(start=LEFT, end=RIGHT)))
gwas_catalog_gr <- with(gwas_catalog, GRanges(seqnames=CHR, ranges=IRanges(start=LEFT, end=RIGHT)))
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

Find overlaps between clumped, COJO, and previous results. Append and
then reduce all regions

``` r
all_gr <- reduce_ranges(bind_ranges(cojo_gr, wray_gr, howard_gr, levey_gr, giannakopoulou_gr, als_gr, meng_gr, gwas_catalog_gr))
```

Find overlaps and make lists for an upset plot. Take the index from the
combined set as the element, then find which of them are found within
each of the sets of SNPs. The upset plot function handles finding each
combination of intersections.

``` r
hits_upset <- list(MDD3_COJO=unique(findOverlaps(all_gr, cojo_gr)@from),
                   Wray=unique(findOverlaps(all_gr, wray_gr)@from),
                   Howard=unique(findOverlaps(all_gr, howard_gr)@from),
                   Levey=unique(findOverlaps(all_gr, levey_gr)@from),
                   Giannakopoulou=unique(findOverlaps(all_gr, giannakopoulou_gr)@from),
                           Als=unique(findOverlaps(all_gr, als_gr)@from),
                   Meng=unique(findOverlaps(all_gr, meng_gr)@from),
                   Catalog=unique(findOverlaps(all_gr, gwas_catalog_gr)@from))
                   
upset(fromList(hits_upset), nsets=8, order.by='freq', text.scale=2)
```

![](/Users/mark/Work/mdd-meta/docs/cojo_files/figure-gfm/upset-1.png)<!-- -->

Find which COJO regions overlap with Howard

``` r
cojo_howard_overlaps <- findOverlaps(cojo_gr, howard_gr)
cojo_howard_overlaps
```

    ## Hits object with 141 hits and 0 metadata columns:
    ##         queryHits subjectHits
    ##         <integer>   <integer>
    ##     [1]         2           1
    ##     [2]         6           2
    ##     [3]         7           3
    ##     [4]        11           4
    ##     [5]        11           5
    ##     ...       ...         ...
    ##   [137]       645          97
    ##   [138]       651          98
    ##   [139]       665         100
    ##   [140]       667         101
    ##   [141]       685         102
    ##   -------
    ##   queryLength: 697 / subjectLength: 102

Count number of regions in Howard that overlap:

``` r
howard %>% slice(unique(cojo_howard_overlaps@to)) %>% count()
```

    ## # A tibble: 1 × 1
    ##       n
    ##   <int>
    ## 1    93

Find which COJO regions overlap with Levey

``` r
cojo_levey_overlaps <- findOverlaps(cojo_gr, levey_gr)
cojo_levey_overlaps
```

    ## Hits object with 299 hits and 0 metadata columns:
    ##         queryHits subjectHits
    ##         <integer>   <integer>
    ##     [1]         2         110
    ##     [2]         6           2
    ##     [3]         7          65
    ##     [4]        10         197
    ##     [5]        11          77
    ##     ...       ...         ...
    ##   [295]       670         154
    ##   [296]       671         154
    ##   [297]       680         120
    ##   [298]       685         153
    ##   [299]       687         211
    ##   -------
    ##   queryLength: 697 / subjectLength: 223

Count number of regions in Levey that overlap:

``` r
levey %>% slice(unique(cojo_levey_overlaps@to)) %>% count()
```

    ## # A tibble: 1 × 1
    ##       n
    ##   <int>
    ## 1   203

Overlaps with previous findings.

``` r
cojo_known_overlaps <- findOverlaps(cojo_gr, reduce_ranges(bind_ranges(wray_gr, howard_gr, levey_gr, giannakopoulou_gr, als_gr, meng_gr, gwas_catalog_gr)))

cojo_known <- cojo %>% slice(unique(cojo_known_overlaps@from))

catalog_known <- rp_gwas_catalog_entries %>% filter(SNP %in% cojo_known$SNP) %>% count(phenotype) %>% arrange(desc(n))
catalog_known
```

    ## # A tibble: 208 × 2
    ##    phenotype                   n
    ##    <chr>                   <int>
    ##  1 Serum_metabolit...        190
    ##  2 Schizophrenia              52
    ##  3 Intelligence_(MTAG)        48
    ##  4 Waist_circumfer...         38
    ##  5 Body_mass_index            34
    ##  6 Depression_(broad)         33
    ##  7 Neuroticism                31
    ##  8 Trans_fatty_acid_levels    30
    ##  9 Autism_spectrum...         29
    ## 10 Glycerophosphol...         27
    ## # ℹ 198 more rows

Newly discovered regions

``` r
cojo_new <- cojo %>% slice(unique(-cojo_known_overlaps@from)) %>% arrange(P)
cojo_new %>%
select(Locus, SNP.Index, CHR, SNP, BP, P, pJ) %>%
group_by(Locus)
```

    ## # A tibble: 298 × 7
    ## # Groups:   Locus [293]
    ##    Locus SNP.Index   CHR SNP               BP        P       pJ
    ##    <dbl>     <dbl> <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ##  1   493         1    14 rs67958427  33778247 3.25e-17 5.70e-17
    ##  2   451         1    12 rs1995514   60789609 6.91e-16 8.85e-16
    ##  3   443         1    12 rs7973992   39168233 4.53e-15 2.92e-14
    ##  4    88         1     2 rs13025467 173797297 8.77e-15 8.14e-15
    ##  5   140         1     3 rs10937241 185822774 1.09e-14 5.47e-15
    ##  6   208         1     5 rs358667   143877756 2.75e-14 2.28e-15
    ##  7   262         1     7 rs10238024   8680443 3.44e-14 6.03e-15
    ##  8   218         1     5 rs3822662  170846151 7.05e-14 4.75e-14
    ##  9   167         1     4 rs13134858 115522306 9.35e-14 2.84e-13
    ## 10   146         1     4 rs556798    17259101 1.31e-13 9.30e-13
    ## # ℹ 288 more rows

``` r
rp_gwas_catalog_entries %>% filter(SNP %in% cojo_new$SNP) %>% count(phenotype) %>% arrange(desc(n)) %>% filter(!phenotype %in% catalog_known$phenotype)
```

    ## # A tibble: 38 × 2
    ##    phenotype                    n
    ##    <chr>                    <int>
    ##  1 Type_2_diabetes             11
    ##  2 Body_mass_index_(adult)      3
    ##  3 Diastolic_blood_pressure     2
    ##  4 Glomerular_filt...           2
    ##  5 Survival_in_col...           2
    ##  6 Adiposity                    1
    ##  7 Aspirin_hydroly...           1
    ##  8 Bladder_cancer               1
    ##  9 Blood_urea_nitr...           1
    ## 10 Body_fat_percentage          1
    ## # ℹ 28 more rows

``` r
rp_genes_dist %>% filter(SNP %in% cojo_new$SNP) %>% group_by(SNP) %>% filter(dist2index == min(dist2index)) %>% ungroup() %>% select(gene) %>% distinct()
```

    ## # A tibble: 262 × 1
    ##    gene        
    ##    <chr>       
    ##  1 NTRK3       
    ##  2 NXPH1       
    ##  3 MAGI2       
    ##  4 LOC101927987
    ##  5 CTNNA2      
    ##  6 LHFPL6      
    ##  7 AQP12A      
    ##  8 FAM19A5     
    ##  9 NTM         
    ## 10 TOX3        
    ## # ℹ 252 more rows

Calculate power for previous versus current GWAS

``` r
# use sum of effective sample sizes and case rate of 50% rather than sum of cases and controls
wray_power <- genpwr.calc(calc='es', model='logistic', ge.interaction=NULL, N=405961, Case.Rate=0.5, k=NULL, MAF=seq(0.005, 0.49, by=0.005), Power=0.8, Alpha=5e-8, True.Model='Additive', Test.Model='Additive')

als_power <- genpwr.calc(calc='es', model='logistic', ge.interaction=NULL, N=938012+100100+24274, Case.Rate=0.5, k=NULL, MAF=seq(0.005, 0.49, by=0.005), Power=0.8, Alpha=5e-8, True.Model='Additive', Test.Model='Additive')

meng_power <- genpwr.calc(calc='es', model='logistic', ge.interaction=NULL, N=700089, Case.Rate=0.5, k=NULL, MAF=seq(0.005, 0.49, by=0.005), Power=0.8, Alpha=5e-8, True.Model='Additive', Test.Model='Additive')

mdd3_power <- genpwr.calc(calc='es', model='logistic', ge.interaction=NULL, N=2*1053442, Case.Rate=0.5, k=NULL, MAF=seq(0.005, 0.49, by=0.005), Power=0.8, Alpha=5e-8, True.Model='Additive', Test.Model='Additive')
```

``` r
frq_u_col <- str_subset(names(cojo), 'FRQ_U')

cojo_known_novel <- bind_rows(
mutate(cojo_known, Association='Known'),
mutate(cojo_new, Association='Novel')) %>%
mutate(BETA=log(OR)) %>%
select(SNP, Association, OR, BETA, FRQ=starts_with('FRQ_U')) %>%
mutate(MAF=if_else(FRQ <= 0.5, true=FRQ, false=1-FRQ))

power_known_novel <- bind_rows(
transmute(als_power, Power='Als', MAF, OR=`OR_at_Alpha_5e-08`),
transmute(meng_power, Power='Meng', MAF, OR=`OR_at_Alpha_5e-08`),
transmute(mdd3_power, Power='MDD3', MAF, OR=`OR_at_Alpha_5e-08`)
)

ggplot(cojo_known_novel, aes(x=MAF, y=exp(abs(BETA)))) + 
geom_point(aes(color=Association)) +
geom_line(mapping=aes(y=OR, linetype=Power), data=power_known_novel, size=1) +
scale_y_continuous('OR', limits=c(1, 1.1))
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

    ## Warning: Removed 3 rows containing missing values (`geom_line()`).

![](/Users/mark/Work/mdd-meta/docs/cojo_files/figure-gfm/cojo_known_novel-1.png)<!-- -->
