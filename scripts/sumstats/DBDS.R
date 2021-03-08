library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(TRUE)

text_gz=args[1]
daner_gz=args[2]
logfile=args[3]

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

sumstats <- read_table2(text_gz,
	col_types=cols(
	  marker = col_character(),
	  Pval = col_double(),
	  Effect = col_double(),
	  rsName = col_character(),
	  MAF_PC = col_double(),
	  eurMAF_PC = col_double(),
	  Info = col_double(),
	  Chrom = col_character(),
	  Pos = col_double(),
	  Amin = col_character(),
	  Amaj = col_character(),
	  TGid = col_character(),
	  HRCid = col_character(),
	  PosB37 = col_double(),
	  TGoa = col_character(),
	  TGea = col_character(),
	  TGeaf = col_double()
	))

Ncases <- 10336
Ncontrols <- 142861

FRQ_A_col <- paste0('FRQ_A_', Ncases)
FRQ_U_col <- paste0('FRQ_U_', Ncontrols)

cat(paste('Formatting daner', '\n'))
daner <-
sumstats %>%
filter(Effect != 1 & Pval != 1) %>%
transmute(CHR=str_replace(Chrom, 'chr', ''), SNP=rsName, BP=PosB37, A1=Amin, A2=Amaj,
		 !!FRQ_A_col:=MAF_PC/100,
		 !!FRQ_U_col:=MAF_PC/100,
		 INFO=Info,
		 OR=Effect,
		 SE=sqrt(log(Effect)^2/qchisq(Pval, df=1, lower.tail=F)),
		 P=Pval)
	

cat(paste('Writing daner', '\n'))   
write_tsv(daner, daner_gz)