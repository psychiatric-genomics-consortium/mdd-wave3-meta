#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting HUNT to daner from $(basename $text_gz)" > $log

 # 1	#CHROM
 # 2	POS
 # 3	ID
 # 4	REF
 # 5	ALT
 # 6	QUAL
 # 7	FILTER 
 # 8	INFO

Rscript scripts/sumstats/HUNT.R $text_gz $daner_gz $log