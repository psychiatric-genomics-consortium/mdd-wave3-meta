# Meta-analysis

The primary meta-analysis workflow has the following steps

1. Pull cohorts' summary statistics in using harmonised filenames.
2. Liftover from HG38 to HG19 when necessary
3. Align summary statistics with a target imputation panel.
4. Setup inputs for Ricopili meta-analysis.
5. Submit Ricopili pipeline.

# Setup

The meta-analyis is run on [LISA](https://geneticcluster.org)

Download [RICOPILI](https://sites.google.com/a/broadinstitute.org/ricopili/), decompress the archive, and move the files to your home directory

```
wget https://sites.google.com/a/broadinstitute.org/ricopili/download/rp_bin.2019_Oct_15.001.tar.gz
tar -xvzf rp_bin.2019_Oct_15.001.tar.gz
mv rp_bin ~/rp_bin
```

Customize RICOPILI installation

```
~/rp_bin/rp_config
```

Customize the file `rp_config.custom.txt` [for LISA](https://docs.google.com/spreadsheets/d/1LhNYIXhFi7yXBC17UkjI1KMzHhKYz0j2hwnJECBGZk4/edit#gid=255132922) and then run 

```
~/rp_bin/rp_config
```

again.

# Interactive node

When running Snakemake rules on LISA, start an interactive node rather than running programs on the login node

```
srun -n 16 -t 1:00:00 --pty bash -il
```

# Adding a new cohort

Cohorts are added to the workflow by creating an entry in under `sumstats` `config.yaml` file with the harmonized key and path to the summary statistics file:

```
sumstats:
	FORMAT_COHORT.POP.hgNN.VERSION: PATH/TO/SUMSTATS.gz
```

where `FORMAT_mdd_COHORT.POP.hgNN.RELEASE` is a key encoding:

- `FORMAT`: format of the summary statistics, one of `daner` or `text`. [daner](https://docs.google.com/document/d/1TWIhr8-qpCXB13WCXcU1_HDio8lC_MeWoAg2jlggrtU/edit) files are left-as-is where as `text` summary statistics will be transformed to daner in another workflow step.
- `COHORT`: short, unique name of the cohort (e.g., `UKBB` for "UK Biobank")
- `POP`: ancestries population for meta-analysis. Should match the ancestries superpopulations used by the Ricopili imputation panel (`afr`, `amr`, `eas`, `eur`, `sas`)
- `hgNN`: Genome build, one of `hg19` or `hg38`.
- `RELEASE`: a release or version specifier for the sumstats from this cohort (usually a sub-cohort, descriptor, or date string).

The `COHORT` and `RELEASE` strings should only contain the characters `[A-Za-z0-9_]` (specifically, the period character `.` should be avoided as it is used to demarcate the different parts of the key).

## Linking a cohorts into the workflow

To link a single cohort: 

```
snakemake -j1 --use-conda resources/sumstats/FORMAT_mdd_COHORT.POP.hgNN.RELEASE.gz
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

that matches the cohort name `COHORT` of the summary stats file (`resources/sumstats/text_mdd_COHORT.POP.hgNN.RELEASE.gz`). The script WILL be executed as

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
echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner

# example: rearrange columns to match daner format
zcat $sumstats | tail -n +2 | awk -v OFS='\t' '{print $1, $3, $2, $4, $5, $8, $9, $10, $13, $14, $15}' >> $daner

# compress daner file
gzip --verbose $daner
```

The workflow will handle the naming of the `INPUT` and `OUTPUT` parameters. Run the rule with

```
snakemake -j1 results/sumstats/daner/daner_mdd_COHORT.POP.hgNN.RELEASE.gz
```

## Liftover

The meta-analysis is conducted against genome assembly [GRCh37](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/) (hg19). Summary statistics on build `hg38` are automatically lifted over to `hg19` as part of the workflow. Assuming that the sumstats file

```
resources/sumstats/FORMAT_mdd_COHORT.POP.hg38.RELEASE.gz
```

exists then it can be converted separately to `hg19` by running

```
snakemake -j1 --use-conda results/sumstats/hg19/daner_mdd_COHORT.POP.hg19.RELEASE.gz
```

## Add cohort to the meta-analysis list

The meta-analysis rule for each ancestries group is in the [`rules/meta.smk`](../rules/meta.smk) Snakemake file.The cohort/release names are in a list at the top of the Snakemake file called `cohorts_POP`, where `POP` is the name ancestries superpopulation group name (lowercase). Add a two element entry to this list of the form `["COHORT"], ["RELEASE"]`. For example, if the current entries are

```
cohorts_eur = [["MDD49", "29w2_20w3_1504"], 
["23andMe", "v7_2_202012"],      
["deCODE", "DEPALL_FINAL_WHEAD"]]
```

and the summary statistics we want to add are for the `GenScot` cohort release version `1215a`, then the updated list would be:

```
cohorts_eur = [["MDD49", "29w2_20w3_1504"], 
["23andMe", "v7_2_202012"],      
["deCODE", "DEPALL_FINAL_WHEAD"],
["GenScot", "1215a"]]
```

# Set-up meta

Run the `dataset_POP` commands to do all the preprocessing of the summary statistics.

```
snakemake -j1 --use-conda dataset_eur
```

Run the QC steps to calculate pairwise LDSC regressions between all the input cohorts
```
snakemake -j16 --use-conda meta_ldsc_sumstats_pairs
```

Then run the rule to generate all of the preprocessing QC information

```
snakemake -j1 --use-conda meta_qc
```

Finally, examine the QC output report: `docs/metaqc.html`.

# Run the meta-analysis

The meta-analysis can be submitted to Ricopili with

```
snakemake -j1 --use-conda results/meta/full_POP_v3.N.M.done 
```

where `POP` is the ancestries group and `v3.N.M` is a version number specifying the number of cohorts included, where `N` is the number of PGC cohorts (analysed from genotype data) and `M` is the number of additional cohorts analysed from summary statitics. For example, the above meta analysis with 29 PGC MDD cohorts and three additional cohorts (23andMe, deCode, and GenScot) would be `v3.29.03`.

Once the full and LOO Ricopili meta-analyses jobs finish, create the final versions of the sumstats for downstream analysis:

```
snakemake -j4 --use-conda postimp_rp_all
```