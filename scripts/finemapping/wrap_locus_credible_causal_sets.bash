### Define infiles

loci=$1
inprefix=$2

while read -r line
do
    echo "Rscript "$PWD"/scripts/finemapping/locus_credible_causal_sets.R ${inprefix}.chr${line}.gz"
done < $loci > $PWD/scripts/finemapping/run_locus_credible_causal_sets.bash

bash $PWD/scripts/finemapping/run_locus_credible_causal_sets.bash

output=$(echo $inprefix | sed 's|finemapping/locus_results|finemapping/locus_credible_causal|g' | sed 's|.rp|.rp.complete|g')
touch $output
