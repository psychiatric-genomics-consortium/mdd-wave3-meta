#!/bin/bash

text_gz=$1
daner_gz=$2

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

echo "CHR SNP BP A1 A2 FRQ_A_23424 FRQ_U_192220 INFO OR SE P" > $daner

zcat $text_gz | tail -n +2 | awk '{print $1, $5, $2, $4, $3, $11, $12, $18, exp($8), $9, $7}' >> $daner

gzip --verbose $daner