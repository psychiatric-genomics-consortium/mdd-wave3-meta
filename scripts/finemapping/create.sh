
inputs=("${@:1:22}")
n=$23
output=$24

## Create finemapping jobs

for input in inputs
do
    chr=$(echo $input | grep -Eo 'chr.[0-9]{1,2}' | sed 's/chr.//g')
    
    python3 resources/finemapping/polyfun/create_finemapper_jobs.py \
        --sumstats $input \
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

for file in ${output}_jobs_*
do
    sed -i 's|finemapper.py|resources/finemapping/polyfun/finemapper.py|g' $file
done

## Run jobs

for chr in {22..1}
do
    sh ${output}_jobs_chr$chr.sh
done
