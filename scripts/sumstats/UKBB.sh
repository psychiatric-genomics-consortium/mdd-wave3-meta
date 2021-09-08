#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting UKBB to daner from $(basename $text_gz)" > $log

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

Nca=$(zcat $text_gz | sed -n '2p' | awk '{{print $6/2}}')
Nco=$(zcat $text_gz | sed -n '2p' | awk '{{print $7/2}}')

echo "Nca: ${Nca}" >> $log
echo "Nco: ${Nco}" >> $log

#  1	#CHROM
#  2	POS
#  3	ID
#  4	A1
#  5	AX
#  6	CASE_ALLELE_CT
#  7	CTRL_ALLELE_CT
#  8	A1_CASE_FREQ
#  9	A1_CTRL_FREQ
# 10	MACH_R2
# 11	FIRTH?
# 12	TEST
# 13	OBS_CT
# 14	OR
# 15	LOG(OR)_SE
# 16	P

echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner

zcat $text_gz | tail -n +2 | awk -v OFS='\t' '{if($1 == "X") {print 23, $3, $2, $4, $5, $8, $9, $10, $14, $15, $16} else {print $1, $3, $2, $4, $5, $8, $9, $10, $14, $15, $16}}' >> $daner

gzip -f --verbose $daner 2>> $log