library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(parallel)

args <- commandArgs(TRUE)

save.image()

stop()

text_gz=args[1]
daner_gz=args[2]
logfile=args[3]

cat(paste('Counting SNPs', '\n'))
snp_count <- length(count.fields(text_gz, comment.char='#'))


Ncases <- max(sumstats$NCAS)
Ncontrols <- max(sumstats$NCON)

FRQ_A_col <- paste0('FRQ_A_', Ncases)
FRQ_U_col <- paste0('FRQ_U_', Ncontrols)

cat(paste('Formatting daner', '\n'))
daner <-
sumstats %>%
mutate(CHR=if_else(`#CHROM` == 'X', true=23, false=as.numeric(`#CHROM`)),
	   OR=exp(BETA)) %>%
select(CHR, SNP=ID, BP=POS, A1=ALT, A2=REF,
		 !!FRQ_A_col:=FCAS, !!FRQ_U_col:=FCON, INFO=R2,
	   OR, SE, P=PVAL, ngt=NGT)

cat(paste('Writing daner', '\n'))   
write_tsv(daner, daner_gz)