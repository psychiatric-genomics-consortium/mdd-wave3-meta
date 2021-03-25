#!/bin/bash

text_gz=$1
daner_gz=$2
log=$3

echo "Filter EXCEED to daner from $(basename $text_gz)" > $log

daner=$(dirname $daner_gz)/$(basename $daner_gz .gz)

zcat $text_gz | awk 'NR ==1 || ($6 >= 0.01 || $7 >= 0.01)' >> $daner

gzip -f --verbose $daner 2>> $log