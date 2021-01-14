#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Converting MVP to daner from $(basename $text_gz)" > $log

# 1	rsid
# 2	CHR:BP
# 3	CHR
# 4	BP
# 5	A1
# 6	A2
# 7	log(OR)
# 8	SE
# 9	P

Rscript scripts/sumstats/MVP.R $text_gz $daner_gz $log