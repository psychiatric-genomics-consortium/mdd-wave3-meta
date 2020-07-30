#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting ALSPAC to daner from $(basename $text_gz)" > $log

    #  1  SNP
    #  2  CHR
    #  3  BP
    #  4  coded_allele
    #  5  noncoded_allele
    #  6  strand_genome
    #  7  BETA (column labelled "OR")
    #  8  SE
    #  9  P
    # 10  n_cases
    # 11  n_controls
    # 12  AF_coded_all
    # 13  AF_coded_cases
    # 14  AF_coded_controls
    # 15  info

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

Nca=$(zcat $text_gz | sed -n '2p' | awk '{{print $10}}')
Nco=$(zcat $text_gz | sed -n '2p' | awk '{{print $11}}')

echo "Nca: ${Nca}" >> $log
echo "Nco: ${Nco}" >> $log

echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner

zcat $text_gz | tail -n +2 | awk -v OFS='\t' '{print $2+0, $1, $3, $4, $5, $13, $14, $15, exp($7), $8, $9}' >> $daner

gzip --verbose $daner 2>> $log