# Analysis Plans

Plans for analyses to include in primary and who is involved:

**Analysts:** Mark Adams, Swapnil Awasthi, Fabian Strait, Oliver Pain, Xueyi Shen, Olga Giannakopoulou, Alex Kwong, Shuyang Yao, Jackson Thorp, Bochao Lin, Gita Pathak, Abigail ter Kuile, Alish Palmos, Karmel Choi, Jonathan Coleman, David Howard, Na Cai, Zac Gerring, Eske Derks   

## Meta-analysis

- GWAS of genotyped cohorts + meta-analysis of same _(Awasthi)_
- Meta-analysis GWAS of genotyped + summary statistics cohorts 
  - Inverse-variance weighted ([Ricopili](https://sites.google.com/a/broadinstitute.org/ricopili/)) _(Adams)_
  - Structural meta analysis based on phenotyping (eHR – self declared – interview based) ([MTAG](https://github.com/omeed-maghzian/mtag)/[GenomicSEM](https://github.com/MichelNivard/GenomicSEM)) _(Adams, Thorp, Cai)_
  - Final primary versions:
    - Full meta (PGC only)
    - noUKBB (PGC only)
    - top 10k (Public)
    - no23am (Public)
    - noUKBB / no23am (Public)
 
## Omics
 
- Positional-mapping / fine-mapping _(Coleman)_
- Gene-based analysis _(Howard)_
- Expression based mapping / eQTL _(Derks, Gerring)_
- Hi-C (chromosome conformation) ([H-MAGMA](https://github.com/thewonlab/H-MAGMA)) _(Pathak)_
- cell type analysis _(Yao)_
- TWAS / gene-set enrichment _(Pain)_
- Gene annotation: _(Coleman, Howard, ter Kuile, Pathak)_
 
## Prediction
 
- SBayesR - estimate joint effect of SNPs / single threshold PRS ([GCTB](https://cnsgenomics.com/software/gctb/#SummaryBayesianAlphabet)) _(Wray)_
- Replication (or possible meta) in other cohorts
- Leave one out PRS – to be provided to each group _(Awasthi)_
- PheWAS of brain and behavior _(Shen, Choi)_
- Mendelian Randomization
- - Two-sample MR of brain and behavior _(Shen, Choi)_
- - latent-causal analysis ([gLVC](https://github.com/lukejoconnor/LCV)) _(Pathak)_
- - MR-PRESSO _(Lin)_
- - [SMR](https://cnsgenomics.com/software/smr/#Overview) (summary-data-based Mendelian Randomization) _(Pathak, Wray)_
- Cross-ancestries prediction _(Giannakopoulou)_
