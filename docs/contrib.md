# Contributing downstream analyses

## Introduction

This project uses the [Snakemake](https://snakemake.readthedocs.io) build system to generate reproducible results. See the [README](../README.md) for installation instructions. [Why Snakemake](https://vincebuffalo.com/blog/2020/03/04/understanding-snakemake.html)?

1. A build system is an explicit way to represent the dependencies between data, code, and results.
2. Compared with [Make](https://www.gnu.org/software/make), Snakemake offers the flexibility of a scripting language and, as an extension of Python, is easy to read and write. 
3. Workflow tools like [bpipe](http://docs.bpipe.org) are more structured for defining how input files map on to output files. Bpipe is good for projects where you have a lot of input files that all need to be processed the same way. With Snakemake, the conceptual focus is more heavily on telling the system what *output* files we want, and the workflow automatically determines the dependencies necessary to create those files. While both systems are very flexible, Snakemake's conceptual scope is a better fit for a project where we are primarily concern are the results for a particular analysis (GWAS of MDD). 
4. [Snakemake integrates with scripts](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts) in the three main open data science languages: [R](https://www.r-project.org), [Python](https://www.r-project.org), and [Julia](https://julialang.org).
5. As an integrated, scalable platform, [Hail](https://hail.is/) are another good choice for reproducible results, but for this project with contributers from many different groups, Snakemake facilitates workflows running on institutional [clusters](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) in addition to [cloud platforms](https://snakemake.readthedocs.io/en/stable/executing/cloud.html). Some of the workflows need to be run on systems were Spark is not available. 

## New analysis workflow

Each workflow is build around a set of [rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html) for generating the output files from all of the input files (resources, data, scripts, and outputs of other rules).

Most downstream analyses will start from the final meta-analysed summary statistics which can be found:

```
results/distribution/daner_pgc_mdd_COHORTS_POP_hg19_v3.NN.MM.gz
```

where 

- `COHORTS` identifies which cohorts are included. The `COHORTS` will most likely be called `full`, meaning that it includes all cohorts that are part of the analysis, or `no23andMe`, meaning all cohorts except for 23andMe.
- `POP` is the ancestries superpopulation (e.g., `eur` for European ancestries)
- `NN` is the number of PGC clinical cohorts included in the meta-analysis
- `MM` is the number of additional cohorts included in the meta-analysis.
- together `v3.NN.MM` is a verion for this analysis that can be used to track each summary statistics file as more cohorts are added to the eventual data freeze. 

For each analysis, create the following three files (replace `analysis` with a meaningful name for your workflow)

- `rules/analysis.smk`: a Snakemake file of rules for your workflow
- `envs/analysis.yaml`: a Conda environment file listing software dependencies that need to be installed to run your rules
- `docs/analysis.md`: a [Markdown](https://guides.github.com/features/mastering-markdown/) with basic documentation about your workflow.

Other important directories that you will use:

- `resources/analysis`: Store any downloaded resources that your analysis relies on here (note, it is possible that some resources will be shared with other analyses, so coordinate with other analysts to determine the best path to store everything)
- `results/analysis`: Directory to store the output of your analysis. Depending on how things are put together, it likely that some of your rules will depend on outputs from other workflows that are stored here. 
- `scripts`: Directory to store scripts used by your rules. Depending on the number of scripts you need, you might make a subdirectory here too.
- `logs/analysis`: A place to store [log files(https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#log-files).

Importantly, the `resources` and `results` directories are excluded from the version control system, so any outputs that need to be version controlled, such as tables, figures, or reports, will need to be stored elsewhere (for example, perhaps under `docs/reports`).

## Your first rule

See the [basic Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html) before continuing to understand the rule syntax

Your first rule will most likely build off of the summary statistics file. Rather than hardcoding the name of the file, we can use [wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards) to generalise the rule and to extract important information from the sumstats filename. We'll call our rule `part1`

```
rule part1:
    input: "results/distribution/daner_pgc_mdd_{cohorts}_{ancestries}_hg19_v{version}.gz"
    output: "results/analysis/part1_{cohorts}_{ancestries}_hg19_v{version}.out"
    shell: "command {input} {output}"
```

If this rule were run with an input called `results/distribution/daner_pgc_mdd_full_eur_hg19_v3.29.08.gz` then during rule execution, the variable `cohorts` would have the value `full`, `ancestries` would have the value `eur`, and `version` would have the value `3.29.08`. The `output` variable would have the value `results/analysis/part1_full_eur_hg19_v3.29.08.out` and the shell command

```
command results/distribution/daner_pgc_mdd_full_eur_hg19_v3.29.08.gz results/analysis/part1_full_eur_hg19_v3.29.08.out
```

would get executed. Instead of running the shell command, we ask Snakemake the execute the rule just by asking for the full path of the requried output file

```
snakemake -j1 results/analysis/part1_full_eur_hg19_v3.29.08.out
```

(where the `-j1` is a required flag to specify the number of precessor cores to use to run the rule).
 
## Linking your workflow into the main Snakefile

Before your first rule is run, your Snakefile needs to be imported by the [`Snakefile`](../Snakefile) in the main directory. Add a line to the main Snakefile:

```
include: "rules/analysis.smk"
```

## Your second rule

A second can be constructed by making the output of the first rule into an input

```
rule part2:
    input: "results/analysis/part1_{cohorts}_{ancestries}_hg19_v{version}.out"
    output: "results/analysis/part1_{cohorts}_{ancestries}_hg19_v{version}.out"
    script: "Rscript ../scripts/analysis.R"
```

This time, the rule calls an [R script](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts). When this happens, an object called `snakemake` will automatically be loaded into the R session and available to use (the path of the input file will be available as the variable `snakemake@input[[1]]`). Note that the path of the script is *relative* to the location of your Snakemake file in `rules/analysis.smk`. 

## Have your rules automatically discover when a new sumstats file is available

## Setting up an environment

Create an [environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually) file using YAML syntax in `envs/analysis.yaml`. You can search for packages available on conda with the command

```
conda search --channel CHANNEL PACKAGE
```

Most bioinformatics packages are available on the `bioconda` channel while others, such as R and CRAN libraries, are the default chanel and the `conda-forge` channel. For example, our `part2` rule requires R, so we might make an environment file `envs/analysis.yaml` with the contents

```
channels:
  - conda-forge
  - bioconda
dependencies:
  - r-base =4.0.2
  - r-essentials
```

to request R version 4.0.2. The rule would be updated to [use the environment](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) as


```
rule part2:
    input: "results/analysis/part1_{cohorts}_{ancestries}_hg19_v{version}.out"
    output: "results/analysis/part1_{cohorts}_{ancestries}_hg19_v{version}.out"
    conda: "../envs/analysis.yaml"
    script: "Rscript ../scripts/analysis.R"
```

(Note again the relative path to `analysis.yaml`).

## Rules to download resources

## Rules to share results
