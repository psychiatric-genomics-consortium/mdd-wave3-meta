#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting AGDS to daner from $(basename $text_gz)" > $log

#  1	SNP
#  2	CHR
#  3	BP
#  4	A1
#  5	A2
#  6	FREQ
#  7	BETA
#  8	SE
#  9	P
# 10	N

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

Nca=12123
Nco=12684

echo "Nca: ${Nca}" >> $log
echo "Nco: ${Nco}" >> $log

echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner

zcat $text_gz | tail -n +2 | awk -v OFS='\t' '{print $2, $1, $3, $4, $5, $6, $6, 1, exp($7), $8, $9}' >> $daner

gzip -f --verbose $daner 2>> $log