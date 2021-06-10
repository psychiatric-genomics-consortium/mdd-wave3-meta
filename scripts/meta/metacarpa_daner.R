library(dplyr)
library(readr)
library(tidyr)

carpa_txt <- snakemake@input[[1]]
carpa <- read_tsv(carpa_txt, col_types=list(`chr:pos`=col_character()))

# find largest sample size
Nco <- Nca <- round(max(carpa$n)/2)

frq_a_col <- paste('FRQ', 'A', Nca, sep='_')
frq_u_col <- paste('FRQ', 'U', Nca, sep='_')

daner_wald <- 
carpa %>%
separate(`chr:pos`, c('CHR', 'BP'), sep=':') %>%
transmute(CHR, SNP=rsid, BP,
A1=effect_allele, A2=neffect_allele,
!!frq_a_col:=effect_allele_frequency, 
!!frq_u_col:=effect_allele_frequency,
INFO=1, OR=exp(beta), SE=se, P=p_wald, Direction=effects, 
Nca=round(n/2), Nco=round(n/2), Neff_half=round(n/2))


daner_gz <- snakemake@output[[1]]
write_tsv(daner_wald, daner_gz)



# lambda_gc <- median(qchisq(P, 1, lower.tail=F)) / qchisq(0.5,1)
