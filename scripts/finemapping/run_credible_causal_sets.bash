### Define infiles

loci=$1
inprefix=$2

while read -r line
do
    echo "Rscript /home/coleman/MDD/mdd-meta/scripts/finemapping/credible_causal_sets.R ${inprefix}.chr${line}.gz"
done < $loci > actually_run_credible_causal_sets.bash

bash actually_run_credible_causal_sets.bash
