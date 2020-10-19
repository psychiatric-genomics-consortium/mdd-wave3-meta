#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting MoBa to daner from $(basename $text_gz)" > $log

#  1	SNP
#  2	A1
#  3	CHR
#  4	BP
#  5	N
#  6	OR
#  7	SE
#  8	P
#  9	A2
# 10	Freq_A
# 11	Freq_U
# 12	INFO

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

Nca=$(zcat $text_gz | head -n 1 | awk '{print $10'} | awk -F_ '{print $3}')
Nco=$(zcat $text_gz | head -n 1 | awk '{print $11'} | awk -F_ '{print $3}')

echo "Nca: ${Nca}" >> $log
echo "Nco: ${Nco}" >> $log

echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner

zcat $text_gz | dos2unix | tail -n +2 | awk -v OFS='\t' '{print $3, $1, $4, $2, $9, $10, $11, $12, $6, $7, $8}' >> $daner

gzip -f --verbose $daner 2>> $log