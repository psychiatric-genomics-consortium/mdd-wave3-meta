
loci=$1
inputprefix=$2
nfile=$(cat $3)
output=$4

## Create finemapping jobs

awk -v output=$output \
    -v inputprefix=$inputprefix \
    -v nfile=$nfile \
    '{print "python3 resources/finemapping/polyfun/finemapper.py --chr "$1" --start "$2" --end "$3" --out "output".chr"$1"."$2"_"$3".gz --method susie --sumstats "inputprefix"."$1".snpvar_ridge_constrained.gz --n "nfile" --memory 1 --max-num-causal 1 --allow-missing"  > output"_jobs_chr"$1".sh"}' $loci 

## Run jobs

for chr in {22..1}
do
    sh ${output}_jobs_chr$chr.sh
done

touch $output
