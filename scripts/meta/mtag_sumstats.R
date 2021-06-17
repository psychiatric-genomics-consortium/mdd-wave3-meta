library(dplyr)
library(readr)
library(stringr)

daner_gz <- snakemake@input[[1]]
daner <- read_table2(daner_gz, col_types=cols("SNP"=col_character()))

# check if sumstats have information on Neff
if("Neff_half" %in% names(daner)) {
	daner_n <- daner %>% mutate(N=2*Neff_half)
} else {
	# get sample sizes from column names
	col_names <- names(daner)
	frq_a_col <- str_subset(col_names, "FRQ_A")
	Nca <- as.numeric(str_extract(frq_a_col, "[[:digit:]]+"))
	frq_u_col <- str_subset(col_names, "FRQ_U")
	Nco <- as.numeric(str_extract(frq_u_col, "[[:digit:]]+"))
	
	# calculate effective sample size
	Neff <- (4*Nca*Nco) / (Nca + Nco)
	daner_n %>% mutate(N=Neff)
}

mtag_sumstats <- daner_n %>%
mutate(BETA=log(OR)) %>%
select(snpid=SNP, chr=CHR, bpos=BP, a1=A1, a2=A2,
       freq=starts_with('FRQ_U'), beta=BETA, se=SE, pval=P, n=N)
	   
sumstats_txt <- snakemake@output[[1]]
write_tsv(mtag_sumstats, sumstats_txt)