Conditional and Joint Analysis
================

``` r
qc <- snakemake@params$qc
```

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

-   Sumstats:
    results/cojo/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.rp.qc.gz  

-   Clump file:
    results/distribution/daner_pgc_mdd_full_eur_hg19_v3.49.24.05.gz.p4.clump.areator.sorted.1mhc  

-   COJO regions: 549  

-   Clumped SNPs: 738  

-   COJO Selected SNPs: 594  

-   Singleton regions: 63  

-   COJO+Clump SNPs: 606  

-   COJO Final SNPs: 543  

-   COJO Final SNPs p \<= 5e-8, pJ \<= 5e-8: 512  

-   COJO Final SNPs p > 5e-8, pJ \<= 5e-8: 29
