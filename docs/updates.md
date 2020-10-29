# Project updates

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
