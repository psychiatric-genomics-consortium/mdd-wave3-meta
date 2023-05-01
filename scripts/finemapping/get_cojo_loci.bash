
input=$1
output_finemapping=$2
output_credible=$3

awk 'NR > 1 {print $4, $25, $26}' $input > $output_finemapping
awk 'NR > 1 {{print $4"."$25"_"$26}}' $input > $output_credible
