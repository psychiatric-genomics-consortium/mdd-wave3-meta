#!/bin/bash

text_gz=$1
daner_gz=$2

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

Nca=$(zcat $text_gz | sed -n '2p' | awk '{{print $6/2}}')
Nco=$(zcat $text_gz | sed -n '2p' | awk '{{print $7/2}}')

echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner

zcat $text_gz | tail -n +2 | awk -v OFS='\t' '{print $1, $3, $2, $4, $5, $8, $9, $10, $13, $14, $15}' >> $daner

gzip --verbose $daner