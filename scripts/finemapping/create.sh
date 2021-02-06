
#python3 -m rpy2.situation

inputprefix=$1
nfile=$2
output=$3

n=$(cat $nfile)

## Create finemapping jobs

for chr in {1..22}
do
    python3 resources/finemapping/polyfun/create_finemapper_jobs.py \
        --sumstats ${inputprefix}.${chr}.snpvar_ridge_constrained.gz \
        --n $n \
        --method susie \
        --max-num-causal 1 \
        --out-prefix $output \
        --jobs-file ${output}_jobs_chr$chr.sh \
        --chr $chr \
        --regions-file resources/finemapping/polyfun/ukb_regions.tsv.gz \
        --allow-missing
done

## Add path to finemapper

# for file in ${output}_jobs_*
# do
#     sed -i 's|finemapper.py|resources/finemapping/polyfun/finemapper.py|g' $file
# done

## Run jobs

for chr in {22..1}
do
    sh ${output}_jobs_chr$chr.sh
done
