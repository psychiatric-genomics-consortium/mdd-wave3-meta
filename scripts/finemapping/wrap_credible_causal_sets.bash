### Define infiles

loci=$1
inprefix=$2

while read -r line
do
    echo "Rscript "$PWD"/scripts/finemapping/credible_causal_sets.R ${inprefix}.chr${line}.gz"
done < $loci > $PWD/scripts/finemapping/run_credible_causal_sets.bash

bash $PWD/scripts/finemapping/run_credible_causal_sets.bash

output=$(echo $inprefix | sed 's|finemapping/results|finemapping/credible_causal|g' | sed 's|.rp|.rp.complete|g')
touch $output