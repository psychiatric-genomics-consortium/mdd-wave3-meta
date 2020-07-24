#!/bin/bash

text_gz=$1
daner_gz=$2

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

Nca=$(zcat $text_gz | sed -n '2p' | awk '{{print $6/2}}')
Nco=$(zcat $text_gz | sed -n '2p' | awk '{{print $7/2}}')

echo "CHR SNP BP A1 A2 FRQ_A_${Nca} FRQ_U_${Nco} INFO OR SE P" > $daner

zcat $text_gz | tail -n +2 | awk '{print $1, $3, $2, $4, $5, $8, $9, $10, $13, $14, $15}' >> $daner

gzip --verbose $daner