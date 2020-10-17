#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting FinnGen to daner from $(basename $text_gz)" > $log

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_23424\tFRQ_U_192220\tINFO\tOR\tSE\tP" > $daner

# convert chromsome X to 23
zcat $text_gz | tail -n +2 | awk -v OFS='\t' '{if($1 == "X") {print 23, $5, $2, $4, $3, $11, $12, $18, exp($8), $9, $7} else {print $1, $5, $2, $4, $3, $11, $12, $18, exp($8), $9, $7}}' >> $daner

gzip -f --verbose $daner 2>> $log