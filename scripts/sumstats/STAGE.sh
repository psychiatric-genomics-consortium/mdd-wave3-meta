#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting SWE STAGE to daner from $(basename $text_gz)" > $log

#  1	CHR
#  2	SNP
#  3	POS
#  4	A1
#  5	A2
#  6	N
#  7	AF1
#  8	BETA
#  9	SE
# 10	P
# 11	INFO

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

# filling in based on preliminary information
Nca=421
Nco=9134

echo "Nca: ${Nca}" >> $log
echo "Nco: ${Nco}" >> $log

echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner

# linear to log-odds transformation
# k = 421 / (421 + 9124) = 0.044
# t = k * (1 - k) = 0.042
# log-odds = beta / t


zcat $text_gz | tail -n +2 | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $7, $7, $11, exp($8/0.042), $9/0.042, $10}' >> $daner

gzip -f --verbose $daner 2>> $log