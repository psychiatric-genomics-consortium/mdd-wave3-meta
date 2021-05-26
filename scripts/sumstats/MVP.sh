#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3
reference_dir=$4

echo "Converting MVP to daner from $(basename $text_gz)" > $log

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)


# release: ICDdep_AllSex_202101
# 1	rsid
# 2	CHR:BP
# 3	CHR
# 4	BP
# 5	A1
# 6	A2
# 7	log(OR)
# 8	SE
# 9	P

if [[ $text_gz =~ ICDdep_AllSex_202101 ]]; then
	Rscript scripts/sumstats/MVP.R $text_gz $daner_gz $log $reference_dir
fi

# release: zeurREL4icd_depFULL
#  1	rsid
#  2	CHROM.POS
#  3	CHROM
#  4	POS
#  5	A1
#  6	A2
#  7	A1_FREQ
#  8	MACH_R2
#  9	TEST
# 10	OBS_CT
# 11	OR
# 12	SE
# 13	P

if [[ $text_gz =~ zeurREL4icd_depFULL ]]; then
    Nca=66993
	Nco=74366
	echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner
	zcat $text_gz | awk -v OFS='\t' 'NR > 1 {print $3, $1, $4, $5, $6, $7, $7, $8, $11, $12, $13}' >> $daner
	gzip -f --verbose $daner 2>> $log
fi