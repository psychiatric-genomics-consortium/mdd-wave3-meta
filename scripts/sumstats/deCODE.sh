#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting deCODE to daner from $(basename $text_gz)" > $log

	#  1	SNP
	#  2	Chr
	#  3	PosB38
	#  4	otherAllele
	#  5	effectAllele
	#  6	eaFrq
	#  7	Beta
	#  8	SE
	#  9	P
	# 10	Info
	# 11	Flip
	# 12	TGid
	# 13	HRCid
	# 14	PosB37
	# 15	TGoa
	# 16	TGea
	# 17	TGeaf
	# 18	MType
	# 19	R2
	# 20	FDist
	# 21	OR

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

# filling in based on preliminary information
Nca=20000
Nco=28000

echo "Nca: ${Nca}" >> $log
echo "Nco: ${Nco}" >> $log

echo -e "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_${Nca}\tFRQ_U_${Nco}\tINFO\tOR\tSE\tP" > $daner

zcat $text_gz | tail -n +2 | awk -v OFS='\t' '{{sub("chr", "", $2)}{print $2, $12, $14, $5, $4, $6, $6, $10, $21, $8, $9}}' >> $daner

gzip -f --verbose $daner 2>> $log