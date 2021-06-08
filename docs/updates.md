# Project updates

# 8 June 2021

- **v3.49.24.02**: Update BioVU to NoCov version with INFO

# 28 May 2021

- **v3.49.24.01**: Fold HUNT cohorts with METACARPA

# 26 May 2021

- **v3.49.24**: Add 2 wave 3 genotype cohorts (MDD47); SHARE; MVP Release 4 _[EUR]_

# 22 Mar 2021

- **v3.47.23** Add 18 wave 3 genotyped cohorts (MDD47) _[EUR]_

# 8 Mar 2021

- **v3.29.23** Add DBDS. _[EUR]_

# 1 Mar 2021

- **v3.29.22** Add Takeda. Update STAGE and 23andMe _[EUR]_

# 2 Feb 2021

- Renamed primary project branch to `main`. Run the following to update the local environment
```
git branch -m master main
git fetch origin
git branch -u origin/main main
```

# 15 Jan 2021

- **v3.29.21** Add BioVU, EXCEED, MVP _[EUR]_
- calculate LDSC rg with MDD29 (PGC clinical cohorts)

# 4 Dec 2020

- **v3.29.18** Add Basic, iCBT. _[EUR]_
- Switch to HRC 1.1 panel for alignment and clumping. 

# 20 Nov 2020

- **v3.29.16** Add AGDS. _[EUR]_

# 19 Nov 2020

- **v3.29.15** Add PREFECT, STAGE. Update iPSYCH, UKBB _[EUR]_

# 6 Nov 2020

- Meta-analyzed sumstats available directly on Lisa. Set-up config file with an entry for `lisa` under the `distribution` settings that contains the path for the distribution files.
  ```
  remote:
    distribution:
      analyst: 
      public: 
      lisa: /path/to/pgcmdd
    ```
  Then run the normal downstream full for the full results
  ```
  snakemake -j1 downstream_full
  ```

# 27 Oct 2020

- **v3.29.13** Add ESTBB, MoBa, Hunt. Update deCODE. _[EUR]_
- **v3.00.02** Add 23andMe, Taiwan. _[EAS]_

# 2 Sep 2020

- **v3.29.10** Add Partners Biobank cohort
- available LOO cohorts for downstream analysis: `noUKBB`, `no23andMe`, `noALSPAC`

# 20 Aug 2020

- `analysis_version` variables defined as another way to refer to the most recent version of the analysis in Snakemake rules.
- `noUKBB` sumstats available for downstream analysis. Run rule
  ```
  snakemake -j1 downstream_noUKBB
  ```

# 19 Aug 2020

- **v3.29.09** Add Airwave cohort.

# 18 Aug 2020

- Development of [project results website](https://psychiatric-genomics-consortium.github.io/mdd-meta/) in the [TWAS branch](https://github.com/psychiatric-genomics-consortium/mdd-meta/tree/twas).
