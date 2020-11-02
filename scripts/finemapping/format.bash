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
<(gunzip -c $input | awk 'NR > 1 {print $1, $2, $3, $4, $5, $7, $8, $9, $10, $11, $12}') \
> $outputroot

# Identify multiallelics for drop

LANG=C fgrep -wf <(awk '{print $2}' $outputroot | sort | uniq -d) \
<(awk '{print $2}' $outputroot) > Duplicate_Names_To_Drop.txt

echo "Duplicate variant names dropped:"
wc -l Duplicate_Names_To_Drop.txt


LANG=C fgrep -wf <(awk '{print $1"_"$3}' $outputroot | sort | uniq -d) \
<(awk '{print $0, $1"_"$3}' $outputroot) | \
awk '{print $2}' > Duplicate_Positions_To_Drop.txt

echo "Duplicate variant positions dropped:"
wc -l Duplicate_Positions_To_Drop.txt

cat Duplicate_Names_To_Drop.txt Duplicate_Positions_To_Drop.txt | sort | uniq > Total_Duplicates_To_Drop.txt

echo "Total duplicates dropped:"
wc -l Total_Duplicates_To_Drop.txt

# Drop any duplicates
LANG=C fgrep -wvf Total_Duplicates_To_Drop.txt $outputroot > $outputroot.noduplicates
mv $outputroot.noduplicates $outputroot 

# Gzip

gzip $outputroot

# Tidy up

rm Duplicate_Names_To_Drop.txt Duplicate_Positions_To_Drop.txt Total_Duplicates_To_Drop.txt $outputroot.noduplicates
