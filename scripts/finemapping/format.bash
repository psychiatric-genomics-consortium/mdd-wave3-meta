### Format sumstats for munge

input=$1
output=$2

outputroot=$(echo $output | sed 's/.gz//g' )

echo "Duplicate variant names:"
gunzip -c $input | awk '{print $2}' | sort | uniq -d | wc -l

echo "Duplicate variant positions:"
gunzip -c $input | awk '{print $1"_"$3}' | sort | uniq -d | wc -l

# Limit to necessary columns, rename FREQ column

cat \
<(echo "CHR SNP BP A1 A2 FREQ INFO OR SE P N") \
<(gunzip -c $input | awk 'NR > 1 {print $1, $2, $3, $4, $5, $7, $8, $9, $10, $11, $17 + $18}') \
> $outputroot

# Identify multiallelics for drop

LANG=C grep -wf <(awk '{print $2}' $outputroot | sort | uniq -d) \
<(awk '{print $2}' $outputroot) > results/finemapping/Duplicate_Names_To_Drop.txt

echo "Duplicate variant names dropped:"
wc -l results/finemapping/Duplicate_Names_To_Drop.txt

LANG=C grep -wf <(awk '{print $1"_"$3}' $outputroot | sort | uniq -d) \
<(awk '{print $2, $1"_"$3}' $outputroot) | \
awk '{print $1}' > results/finemapping/Duplicate_Positions_To_Drop.txt

echo "Duplicate variant positions dropped:"
wc -l results/finemapping/Duplicate_Positions_To_Drop.txt

cat results/finemapping/Duplicate_Names_To_Drop.txt results/finemapping/Duplicate_Positions_To_Drop.txt | sort | uniq > results/finemapping/Total_Duplicates_To_Drop.txt

echo "Total duplicates dropped:"
wc -l results/finemapping/Total_Duplicates_To_Drop.txt

# Drop any duplicates
LANG=C fgrep -wvf results/finemapping/Total_Duplicates_To_Drop.txt $outputroot > $outputroot.noduplicates
mv $outputroot.noduplicates $outputroot 

# Gzip

gzip $outputroot

