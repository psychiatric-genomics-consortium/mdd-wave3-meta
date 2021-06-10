#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting SWE STAGE to daner from $(basename $text_gz)" > $log

#  1	CHR: chromosome
#  2	POS: genome position 
#  3	SNPID: variant ID
#  4	Allele1: allele 1
#  5	Allele2: allele 2
#  6	AC_Allele2: allele count of allele 2
#  7	AF_Allele2: allele frequency of allele 2
#  8	N: sample size
#  9	BETA: effect size of allele 2
# 10	SE: standard error of BETA
# 11	Tstat: score statistic of allele 2
# 12	p.value: p value (with saddlepoint approximation(SPA) applied for binary traits)
# 13	p.value.NA: p value when SPA is not applied 
# 14	Is.SPA.converge: whether SPA is converged or not 
# 15	varT: estimated variance of score statistic with sample relatedness incorporated
# 16	varTstar: variance of score statistic without sample relatedness incorporated
# 17	AF.Cases: allele frequency of allele 2 in cases 
# 18	AF.Controls: allele frequency of allele 2 in controls 
# 19	INFO: info score (only matched with INFO>0.1)


daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

# sample sizes

Nca=421
Nco=9134

echo "Nca: ${Nca}" >> $log
echo "Nco: ${Nco}" >> $log

echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner

zcat $text_gz | tail -n +2 | awk -v OFS='\t' '{if($19 != "NA"){print $1, $3, $2, $5, $4, $7, $7, $19, exp($9), $10, $12}}' >> $daner


gzip -f --verbose $daner 2>> $log
