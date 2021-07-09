# Conditional and Joint Analysis

We ran a [conditional and joint analysis](https://www.nature.com/articles/ng.2213) using [GCTA](https://cnsgenomics.com/software/gcta/#COJO) to refine the list of independent loci.

- Final meta-analysed SNPs from the Ricopili pipeline were used.
- Ricopili was used for initial clumping with index SNPs identified with $p < 0.0001$ and $r^2 < 0.1$ within 3000kb windows. The extended MHC region was clumped as a single region.
- Sumstats were filtered for MAF >= 0.01 and INFO > 0.6
- Regions with a genome-wide significant SNP (p < 5-e8) were identified from the clumped results. Regions within 50kb of each other were merged.
- SNPs from these regions were extracted filtered to unrelated of self- and genotype-identified European ancestry participants from UK Biobank.
- A conditional analysis was performed on each region using the filtered sumstats superimposed on the UK Biobank LD structure.
