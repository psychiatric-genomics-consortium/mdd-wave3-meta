
input=$1
output_finemapping=$2
output_credible=$3

awk 'NR > 1 {print $3, $26, $27}' $input > $output_finemapping
awk 'NR > 1 {{print $3"."$26"_"$27}}' $input > $output_credible
