# Meta-analysis

The primary meta-analysis workflow has the following steps

1. Pull cohorts' summary statistics in using harmonised filenames.
2. Liftover from HG38 to HG19 when necessary
3. Align summary statistics with a target imputation panel.
4. Setup inputs for Ricopili meta-analysis.
5. Submit Ricopili pipeline.

## Adding a new cohort

Cohorts are added to the workflow by creating an entry in under `sumstats` `config.yaml` file with the harmonized key and path to the summary statistics file:

```
sumstats:
	FORMAT_COHORT.POP.hgNN.VERSION: PATH/TO/SUMSTATS.gz
```

where `FORMAT_mdd_COHORT.POP.hgNN.VERSION` is a key encoding:

- `FORMAT`: format of the summary statistics, one of `daner` or `text`. [daner](https://docs.google.com/document/d/1TWIhr8-qpCXB13WCXcU1_HDio8lC_MeWoAg2jlggrtU/edit) files are left-as-is where as `text` summary statistics will be transformed to daner in another workflow step.
- `COHORT`: short, unique name of the cohort (e.g., `UKBB` for "UK Biobank")
- `POP`: ancestries population for meta-analysis. Should match the ancestries superpopulations used by the Ricopili imputation panel (`afr`, `amr`, `eas`, `eur`, `sas`)
- `hgNN`: Genome build, one of `hg19` or `hg38`.
- `VERSION`: a version specifier for the sumstats from this cohort (usually a version number or date string).

The `COHORT` and `VERSION` strings should only contain the characters `[A-Za-z0-9_]` (specifically, the period character `.` should be avoided as it is used to demarcate the different parts of the key).

## Linking a cohorts into the workflow

To link a single cohort: 

```
snakemake -j1 resources/sumstats/FORMAT_mdd_COHORT.POP.hgNN.VERSION.gz
```

To link all cohort summary statistics listed in `config.yaml`:

```
snakemake -j1 sumstats
```

## Formatting sumstats to daner

For summary statistics that are not already in daner format, create a script called

```
scripts/sumstats/COHORT.sh
```

that matches the cohort name `COHORT` of the summary stats file (`resources/sumstats/text_mdd_COHORT.POP.hgNN.VERSION.gz`). The script WILL be executed as

```
sh scripts/sumstats/COHORT.sh INPUT OUTPUT
```

where the script takes the first argument `INPUT` as the name of the gzip-compressed summary statistics text file and write out a gzip-compressed, daner-formatted file to `OUTPUT`. A general pattern for such scripts is:

```
#!/bin/bash

# first argument is text.gz sumstats input file
text_gz=$1
# second argument is daner.gz output file
daner_gz=$2

# location of output file without the .gz extension
daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

# example: determine number of cases and controls from columns 6 and 7
Nca=$(zcat $sumstats | sed -n '2p' | awk '{{print $6/2}}')
Nco=$(zcat $sumstats | sed -n '2p' | awk '{{print $7/2}}')

# write out daner header line (tab-separated)
echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_23424\tFRQ_U_192220\tINFO\tOR\tSE\tP" > $daner

# example: rearrange columns to match daner format
zcat $sumstats | tail -n +2 | awk -v OFS='\t' '{print $1, $3, $2, $4, $5, $8, $9, $10, $13, $14, $15}' >> $daner

# compress daner file
gzip --verbose $daner
```

The workflow will handle the naming of the `INPUT` and `OUTPUT` parameters. Run the rule with

```
snakemake -j1 results/sumstats/daner/daner_mdd_COHORT.POP.hgNN.VERSION.gz
```