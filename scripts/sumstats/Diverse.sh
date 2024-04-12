#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3
reference_dir=$4

echo "Converting Diverse MA to daner from $(basename $text_gz)" > $log

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

Rscript scripts/sumstats/Diverse.R $text_gz $daner_gz $log $reference_dir
