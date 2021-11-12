input=$1
output=$2


declare -i cases controls n
cases=$(gunzip -c $input | head -1 | awk '{{print $6}}' | sed 's/FRQ_A_//g')
controls=$(gunzip -c $input | head -1 | awk '{{print $7}}' | sed 's/FRQ_U_//g')
n=$cases+$controls
echo $n > $output
