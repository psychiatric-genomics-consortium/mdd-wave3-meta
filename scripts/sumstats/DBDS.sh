#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting Danish Blood Donor Study to daner from $(basename $text_gz)" > $log

#  1	marker
#  2	Pval
#  3	Effect (as OR)
#  4	rsName
#  5	MAF_PC (MAF %)
#  6	eurMAF_PC
#  7	Info
#  8	Chrom (as "chrN")
#  9	Pos
# 10	Amin
# 11	Amaj
# 12	TGid
# 13	HRCid
# 14	PosB37
# 15	TGoa
# 16	TGea
# 17	TGeaf

Rscript scripts/sumstats/DBDS.R $text_gz $daner_gz $log
