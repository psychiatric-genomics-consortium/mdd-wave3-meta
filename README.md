# PGC MDD3 Meta-analysis

Working toward the next meta-analysis ("MDD3") by the Major Depressive Disorder Working Group of the Psychiatric Genomics Consortium.

## Embargo date

These data are private to MDD Working Group. All results found here cannot be share, discussed, or presented in any way without explicit permission from the Working Group chairs. 

## Project overview

Meta-analysis of cohorts from the [MDD2](https://doi.org/10.1038/s41588-018-0090-3) and [follow-up and replication](https://doi.org/10.1038/s41593-018-0326-7) genome-wide association studies of Major Depressive Disorder:

- PGC clinical cohorts ("MDD29")
- [deCODE](http://www.decode.com)
- [Generation Scotland](https://www.ed.ac.uk/generation-scotland/)
- [GERA](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000674.v1.p1)
- [iPsych](https://ipsych.dk)
- [UK Biobank](https://www.ukbiobank.ac.uk)
- [23andMe](https://www.23andme.com/) (v5)

and the following new cohorts:

- [ALSPAC](http://www.bristol.ac.uk/alspac/)
- [FinnGen](https://www.finngen.fi/en)
- [23andMe](https://www.23andme.com/) (v7.2)

…plus more as they are incorporated into the analysis.

Analysis conducted on [LISA](https://geneticcluster.org).

### Step 1

Clone the repository

```
git clone git@github.com:psychiatric-genomics-consortium/mdd-meta.git
cd mdd-meta
```


### Step 2


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

### Step 3

Install [Anaconda](https://www.anaconda.com).

```
srun -n 16 -t 1:00:00 --pty bash -il
wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
sh Anaconda3-2020.02-Linux-x86_64.sh
```

Install dependencies and load the environment

```
conda env create --file environment.yaml
conda activate mdd3
```

### Step 4

Configure the analysis workflow. Make a copy of the configuration file

```
cp config.yaml-template config.yaml
```

Then edit and fill in ```config.yaml``. Under the `daner` entry, fill in the `PATH` to each daner file.

### Step 4

Run the meta analysis

```
```


## Checking the results

## Built With

- [RICOPILI](https://sites.google.com/a/broadinstitute.org/ricopili)

## Analysts

* **Mark James Adams** - *analyst* - [Edinburgh](https://www.ed.ac.uk)
* **Swapnil Awasthi** - *analyst* - [Broad](https://www.broadinstitute.org/)
* **Fabian Strait** - *analyst* - [CIMH](https://www.zi-mannheim.de/)
* **Olga Giannakopoulou** - *analyst* [UCL](http://www.bristol.ac.uk/alspac/)
* **David Howard** - *analyst* - [KCL](https://www.kcl.ac.uk/)
* **Karoline Kuchenbaecker** - *group lead* [UCL](http://www.bristol.ac.uk/alspac/)
* **Naomi Wray** - *analytical group director* - [Queensland](https://cnsgenomics.com/)
* **Stephan Ripke** - *analytical group director* - [Broad](https://www.broadinstitute.org/)
* **Cathryn Lewis** - *workgroup chair* - [KCL](https://www.kcl.ac.uk/)
* **Andrew McIntosh** - *workgroup chair* - [Edinburgh](https://www.ed.ac.uk)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE) file for details

## Acknowledgments


