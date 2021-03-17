#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting BioVU to daner from $(basename $text_gz)" > $log

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

# 1	    CHR
# 2	    POS
# 3	    SNPID
# 4	    Allele1
# 5	    Allele2
# 6	    AC_Allele2
# 7	    AF_Allele2
# 8	    imputationInfo
# 9	    N
# 10	BETA
# 11	SE
# 12	Tstat
# 13	p.value
# 14	p.value.NA
# 15	Is.SPA.converge
# 16	varT
# 17	varTstar
# 18	AF.Cases
# 19	AF.Controls

Nca=7757
Nco=24723

echo "Nca: ${Nca}" >> $log
echo "Nco: ${Nco}" >> $log

echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner

zcat $text_gz | awk -v OFS='\t' '{if($1 != "CHR"){print $1, $3, $2, $5, $4, $18, $19, 1, exp($10), $11, $13}}' >> $daner

gzip -f --verbose $daner 2>> $log