# Contributing downstream analyses

## Introduction

This project uses the [Snakemake](https://snakemake.readthedocs.io) build system to generate reproducible results. See the [README](../README.md) for installation instructions. [Why Snakemake](https://vincebuffalo.com/blog/2020/03/04/understanding-snakemake.html)?

1. A build system is an explicit way to represent the dependencies between data, code, and results.
2. Compared with [Make](https://www.gnu.org/software/make), Snakemake offers the flexibility of a scripting language and, as an extension of Python, is easy to read and write. 
3. Workflow tools like [bpipe](http://docs.bpipe.org) are more structured for defining how input files map on to output files. Bpipe is good for projects where you have a lot of input files that all need to be processed the same way. With Snakemake, the conceptual focus is more heavily on telling the system what *output* files we want, and the workflow automatically determines the dependencies necessary to create those files. While both systems are very flexible, Snakemake's conceptual scope is a better fit for a project where we are primarily concerned with the results for a particular analysis (GWAS of MDD). 
4. [Snakemake integrates with scripts](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts) in the three main open data science languages: [R](https://www.r-project.org), [Python](https://www.r-project.org), and [Julia](https://julialang.org).

## Version control

Start by making a new [branch](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-branches) for your `analysis`. This makes it easier to make changes across the whole project without interfering in other work. 

```
git branch analysis
git checkout analysis
```

Keep your branch up-to-date with the main branch:

```
git pull
git push origin analysis
git push --set-upstream origin analysis
git merge main
```

Finally, when you are ready to merge your changes back into the main branch

```
git checkout main
git merge analysis
```

## New analysis workflow

Each workflow is build around a set of [rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html) for generating the output files from all of the input files (resources, data, scripts, and outputs of other rules).

Most downstream analyses will start from the final meta-analysed summary statistics which can be found:

```
results/distribution/daner_pgc_mdd_COHORTS_POP_hg19_v3.NN.MM.RR.rp.gz
```

where 

- `COHORTS` identifies which cohorts are included. The `COHORTS` will most likely be called `full`, meaning that it includes all cohorts that are part of the analysis, or `no23andMe`, meaning all cohorts except for 23andMe.
- `POP` is the ancestries superpopulation (e.g., `eur` for European ancestries)
- `NN` is the number of PGC clinical cohorts included in the meta-analysis
- `MM` is the number of additional cohorts included in the meta-analysis.
- `RR` is the minor revision number
- together `v3.NN.MM.RR` is a version for this analysis that can be used to track each summary statistics file as more cohorts are added to the eventual data freeze. 

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

## Planning out your analysis workflow

Each step in your analysis should be broken down into a "[rule](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)" that maps one or more input files on to one or more output files. Each rule represents a discrete step in your analysis and should be based around running a single command line program or script.

## Your first rule

See the [basic Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html) before continuing to understand the rule syntax

Your first rule will most likely build off of the summary statistics file. Rather than hardcoding the name of the file, we can use [wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards) to generalise the rule and to extract important information from the sumstats filename. We'll call our rule `analysis_part1`

```
rule analysis_part1:
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
rule analysis_part2:
    input: "results/analysis/part1_{cohorts}_{ancestries}_hg19_v{version}.out"
    output: "results/analysis/part2_{cohorts}_{ancestries}_hg19_v{version}.out"
    script: "../scripts/analysis.R"
```

This time, the rule calls an [R script](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts). When this happens, an object called `snakemake` will automatically be loaded into the R session and available to use (the path of the input file will be available as the variable `snakemake@input[[1]]`). Note that the path of the script is *relative* to the location of your Snakemake file in `rules/analysis.smk`. 

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
rule analysis_part2:
    input: "results/analysis/part1_{cohorts}_{ancestries}_hg19_v{version}.out"
    output: "results/analysis/part2_{cohorts}_{ancestries}_hg19_v{version}.out"
    conda: "../envs/analysis.yaml"
    script: "../scripts/analysis.R"
```

(Note again the relative path to `analysis.yaml`). Once the environment is set up, it can be used with the rule by adding the `--use-conda` flag when running Snakemake

```
snakemake -j1 --use-conda results/analysis/part2_full_eur_hg19_v3.29.08.out
```

To test your rule's environment during development:

```
conda env create --name analysis --file envs/analysis.yaml
conda activate analysis
```

To update the environment during development

```
conda env update --name analysis --file envs/analysis.yaml
```

## Have your rules automatically run without specifying the exact name of the output file

The [`expand()`](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-expand-function) function combined with global substitution can be used to automatically run rules on matched set of files. First, `glob_wildcards()` is used to find every file that matches a downloaded sumstats file with the pattern `results/distribution/daner_pgc_mdd_*.gz`, with the matched pattern being stored in the wildcard `sumstats`. We then add another rule, which we simply call `analysis`, that will `expand()` the matched pattern to ask for all files called `results/analysis/part2_*.out`. This will automatically match the output pattern of rule `part2`, and will result in `part1` being run followed by `part2`. 

```
analysis_sumstats_gz, = glob_wildcards("results/distribution/daner_pgc_mdd_{sumstats}.gz")

rule analysis:
    input: expand("results/analysis/part2_{sumstats}.out", sumstats=analysis_sumstats_gz)
```

Once these rules are set up, the whole workflow can by run just by asking for `analysis` by name:

```
snakemake -j1 --use-conda analysis
```

It's also possible to write a rule that will run on only the most recent version of the summary statistics, rather than all of the versions that have been downloaded. This relies on the `analysis_version` variable defined in the `rules/meta.smk` file. `analysis_version` is a list with the value `["3.N.M.R"]` and can be used in `expand()` statements:

```
rule analysis:
  input: expand("results/analysis/part2_full_eur_hg19_v{version}.out", version=analysis_version)
```


## Rules to download resources

Use the [HTTP(S) remote provider](https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html#read-only-web-http-s). The main `Snakemake` file already creates a provider called `HTTP`:

```
rule analysis_download:
  input: HTTP.remote("https://example.com/file.txt")
  output: "resources/analysis/file.txt"
  shell: "cp {input} {output}"
```

## Rules to share results

# Analyses based on cohort-level summary statistics

If your analysis requires individual cohort summary statistics (as opposed to the final meta-analysis summary statistics) then the analysis must be conducted on [LISA](https://geneticcluster.org). 

# Running workflows on a cluster

Snakemake has features for [cluster execution](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). On LISA, workflow stages can be submitted to the batch system with

```
snakemake -jNN --use-conda --cluster 'sbatch -t MM' OUTPUT_FILE_NAME
```

where _`NN`_ is the number of stages that will be submitted to the queue in parallel and _`MM`_ is the runtime allocation for each stage. Other flags can be passed to the cluster as part of the `sbatch` command. When a job is submitted to LISA, an entire node with 16 cores is dedicated to the task. [Job groupings](https://snakemake.readthedocs.io/en/stable/executing/grouping.html) should be used to make 16 be submitted to each node to run in parallel:

```
snakemake -j4 --use-conda --cluster 'sbatch -t 60 -n 16' analysis --groups analysis_part1=group0 analysis_part2=group0 analysis=group1 --group-components group0=16 group1=16
```

will submit jobs in groups of `16` to the cluster while utilising up to `4` nodes running at the same time.

## Running the Snakemake process itself as a job

On LISA, the Snakemake command can be wrapped and passed to the batch system

```
sbatch -t 360 -n 1 --wrap  "snakemake -j4 --use-conda --cluster 'sbatch -t 60 -n 16' OUTFILE"
```

With this command the Snakemake process will be allowed to for up to 6 hours while each stage of the workflow will be allowed 60 minutes. 

## Grouping stages together

By default the `--cluster` command will submit each stage of the workflow as its own job to the batch system. Stages can be grouped together using the `--groups` and `--group-component` flags at the end of the Snakemake command. For example, if `stage1` has many pieces to it, 16 can be run in parallel a single LISA node with

```
sbatch -t 360 -n 1 --wrap "snakemake -j4 --use-conda --cluster 'sbatch -t 60 -n 16' FILENAME.out --groups stage1=group0 --group-components group0=16"
```

When deciding which jobs to group together, it can be useful to get a list of which ones will be executed using a dry run (`-n` flag)

```
snakemake -j1 -n FILENAME.out
```

### Example

For example, grouping jobs together for the meta-analysis pre-processing. The staging, conversion, and alignment stages are grouped together and run as batches of 5 (to avoid using up too much memory) while the LDSC stages can be run in batches of 16.

```
sbatch -t 360 -n 1 --wrap "snakemake -j4 --use-conda --cluster 'sbatch -t 120 -n 16' results/meta/dataset_full_eur_v3.47.23 --groups stage_sumstats=group0 text2daner=group0 hg19=group0 daner=group0 align=group0 meta=group1 meta_ldsc_mdd2=group1 meta_ldsc_munge=group1 dataset_eur=group2 --group-components group0=5 group1=16 group2=1"
```

## Using scratch space

If a job needs to write intermediate files before producing the final output, use "[shadow rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#shadow-rules)"

```
sbatch --wrap "snakemake --shadow-prefix /scratch --cluster 'sbatch -t 60 -n 16' FILENAME.out"
```


## Job scripts

One common issue on clusters (not LISA) for running Snakemake is that the conda environment is not activated on each worker node. In this case it is necessary to make a custom job script that will setup conda. Create a file such as `resources/jobscript.sh` with the contents like:

```
#!/bin/sh
# properties = {properties}

source ~/.bashrc
conda activate base

{exec_job}
```

Then invoke it using the `--jobscript` flag

```
snakemake -j32 --use-conda --cluster 'qsub' --jobscript resources/jobscript.sh OUTPUT_FILE_NAME
```
