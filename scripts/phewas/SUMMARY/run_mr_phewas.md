MR phewas pipeline
================
X Shen
05 January, 2022

-----

### Prepare sumstats

Traits that were significantly associated with MDD PRS are included in
this analysis. PRS PheWAS should be finished before the MR PheWAS.

#### A list of traits for MR PheWAS

Inputs:

  - Any PRS PheWAS output (.rds)
  - Any data dictionary of the final tested fields (.txt)

Output:

  - data/MR/MR\_pheno\_ls\_UKB.RData

<!-- end list -->

``` bash
Rscript scripts/phewas/PREP.MR/PREP.pheno_of_interest.R \
results/phewas/phewas_out_Body_MRI.rds \
results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt \
data/MR/MR_pheno_ls_UKB.RData
```

#### Download sumstats

Inputs:

  - List of phenotypes for MR produced by the previous step
  - Neale’s sumstats spreadsheet
  - BIG 40 sumstats spreadsheet
  - A file for N count for each phenotype in UKB (created by the PRS
    PheWAS pipeline)

Output:

  - Sumstats downloaded in data/MR/MR\_sumstats

If folders weren’t setup yet, run scripts below:

``` bash
mkdir $myscratch/MR_sumstats
ln -s $myscratch/MR_sumstats data/MR/MR_sumstats
mkdir data/MR/MR_sumstats/meta
```

Then run the main script:

``` bash
Rscript scripts/phewas/PREP.MR/PREP.gwas_sumstats_download.R \
data/MR/MR_pheno_ls_UKB.RData \
data/MR/UKBB_GWAS_Imputed_v3_201807.tsv \
data/MR/BIG40_GWAS.csv \
data/phe_count.rds \
data/MR/MR_sumstats
```

It takes roughly 6 hours to run. Consider submitting the job script:

``` bash
qsub scripts/phewas/PREP.MR/job.download_sumstats.sh
```

#### Meta-analyse left & right separated phentypes

Prepare sumstats for meta-analysis:

``` bash
qsub scripts/phewas/PREP.MR/job.reformatgwas.sh
```

The following step needs to be run on an interactive node.

To allow interactive node run without keeping the window open and
retrieve the closed windows, use ‘screen’.

``` bash
screen
cd $pgc3
# Set the mdd_meta folder as the WD
```

Run meta analysis using Mtag. Open multiple screens to run the following
commands to reduce running time.

``` bash
#module load module load igmm/apps/R/4.1.0
#module load anaconda/5.3.1
#conda config --add envs_dirs /exports/igmm/eddie/GenScotDepression/shen/anaconda3/envs/
#conda config --add pkgs_dirs /exports/igmm/eddie/GenScotDepression/shen/anaconda3/pkgs

conda activate mtag

touch data/MR/running_mtag

while read -r a b c || [ -n "$a" ]; do
  if (( $(bc <<< "$(grep $c data/MR/running_mtag | wc -l)==0")  )) && [ ! -f "${c}_mtag_meta.txt" ]; then
    echo $c >> data/MR/running_mtag

    fname="$(basename $c)"
    echo $fname
    bash scripts/phewas/PREP.MR/func.mtag_meta.sh ${a} ${b} ${c} 
  fi
done < data/MR/mtag_list
```

### Local GWAS

#### Set up folders

``` bash
mkdir $myscratch/phewas_gwas
```

#### Prepare phenotypes and GWAS covariates

Inputs:

  - Inputs 1-6: phenotype data

  - Additional covariates

  - A data dictionary file
    (e.g. results/phewas/data\_dictionary/fields.final.brain\_imaging\_QC\_cov\_phenotype.txt)

  - PRS

  - List of phenotypes for local gwas

Output:

  - Directory for outputs: /exports/eddie/scratch/xshen33/phewas\_gwas

<!-- end list -->

``` bash
Rscript scripts/phewas/PREP.MR/PREP.ukb_gwas_pheno_cov.R \
data/dat.imaging_chunk.rds \
data/dat.cognition_chunk.rds \
data/dat.diet_chunk.rds \
data/dat.activity_chunk.rds \
data/dat.mental_health_chunk.rds \
data/dat.loose_field_chunk.rds \
data/dat.addional_covariates_chunk.rds \
results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt \
data/PRS_all.rds \
data/MR/local.gwas.rds \
/exports/eddie/scratch/xshen33/phewas_gwas
```

#### Create pgen format UKB genetic data

This step may be omitted if it’s already been created. QC criteria:

  - MAF \> 0.001

  - MAC \> 20

  - Max alleles == 2

  - Duplicate variants removed

<!-- end list -->

``` bash
qsub scripts/phewas/PREP.MR/job.pgen_format_geneticDat.sh
```

Merge genetic data (if
necessary)

``` bash
qsub scripts/phewas/PREP.MR/job.merge_geneticDat.sh
```

#### Run regenie scripts

``` bash
# regenie Step 1 ==============================================================
while read -r a b c d || [ -n "$a" ]; do

jobname="$(echo $(basename ${a}) | sed "s/pheno_//g" | sed "s/.tsv//g")_step1"
fname="${d}_pred.ls"

if [ ! -f ${fname} ]
   then
     qsub -N $jobname scripts/phewas/PREP.MR/regenie_jobs/func.regenie_step1_phewas_gwas.sh "${a}" "${b}" "${c}" "${d}"
     echo $jobname
fi

done < data/MR/regenie_step1_ls.tsv


# regenie Step 2 ==============================================================

while read -r a b c d || [ -n "$a" ]; do
   jobname="$(echo $(basename ${d}) | sed "s/pheno_//g" )"
   prevjob="$(echo $(basename ${d}) | sed "s/pheno_//g"  | sed "s/_step2/_step1/g")"
   echo $jobname
   echo $prevjob

   qsub -N "${jobname}" -hold_jid "${prevjob}" scripts/phewas/PREP.MR/regenie_jobs/func.regenie_step2_phewas_gwas.sh $a $b $c $d

done < data/MR/regenie_step2_ls.tsv
```

Some GWAS may not finish running before the node exits automatically
after running for 48 hrs. For these ones (check
scripts/phewas/PREP.MR/run\_regenie.sh):

1.  Update a new list of unfinished GWAS to run with their according
    starting block
    
    ``` bash
    Rscript scripts/phewas/PREP.MR/PREP.ukb_gwas_pheno_cov.R \
    data/dat.imaging_chunk.rds \
    data/dat.cognition_chunk.rds \
    data/dat.diet_chunk.rds \
    data/dat.activity_chunk.rds \
    data/dat.mental_health_chunk.rds \
    data/dat.loose_field_chunk.rds \
    data/dat.addional_covariates_chunk.rds \
    results/phewas/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt \
    data/PRS_all.rds \
    data/MR/local.gwas.rds \
    /exports/eddie/scratch/xshen33/phewas_gwas
    ```

2.  Run regenie
    step2
    
    ``` bash
    # regenie Step 2 (remaining) ==================================================
    
    while read -r a b c d f || [ -n "$a" ]; do
       jobname="$(echo $(basename ${d}) | sed "s/pheno_//g" )"
       echo $jobname
       echo "starting block:" $f
    
       qsub -N "${jobname}" scripts/phewas/PREP.MR/regenie_jobs/func.regenie_step2_phewas_gwas_remaining.sh $a $b $c $d $f
    
    done < data/MR/regenie_step2_ls_remaining.tsv
    ```

#### Reformat regenie sumstats

``` bash
qsub scripts/phewas/PREP.MR/job.reformat_regenie.sh
```

-----

### Run MR main analysis

#### Prepare MDD sumstats as exposure

``` bash
bash scripts/phewas/ANALY.MR/job1.prep_mdd_exposure.sh
```

#### Prepare other phenotypes as exposure and outcome

Submit this job script multiple times to reduce running time. It takes
about 1 hr when 30 jobs are running simultaneously.

``` bash
qsub scripts/phewas/ANALY.MR/job2.prep_gwas_expo_outc.sh
```

#### Prepare mdd as outcome

In this step, the script uses the outputs from the previous step. Submit
this job script multiple times to reduce running time.

``` bash
qsub scripts/phewas/ANALY.MR/job3.prep_mdd_outcome.sh
```

#### Create a list of MR analyses to run

``` bash
qsub scripts/phewas/ANALY.MR/job4.mr_input_list.sh
```

#### Run MR analysis

Takes about half an hour.

``` bash
qsub scripts/phewas/ANALY.MR/job5.mr_ukb.sh
```

-----

### Run MR Base

``` bash
qsub scripts/phewas/ANALY.MR/job5.mr_base.sh
```

-----

### Run top SNP MR (sensitivity analysis)

Run from job1 to job5 in folder:
***scripts/phewas/ANALY.MR/ANALY.top\_instrument\_MR***.

Steps are similar with main analysis.

-----

### Run MR using the Howard et al MDD sumstats

Run from job1 to job5 in folder:
***scripts/phewas/ANALY.MR/ANALY.Howard\_etal\_MR***.

Steps are similar with main analysis.
