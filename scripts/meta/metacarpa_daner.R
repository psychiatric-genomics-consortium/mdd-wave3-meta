library(dplyr)
library(readr)
library(tidyr)

mc <- read_tsv('results/meta/metacarpa/pgc_mdd_full_eur_hg19_v3.49.24.mc.txt', col_types=list(`chr:pos`=col_character()))

daner_wald <- 
mc %>%
separate(`chr:pos`, c('CHR', 'BP')) %>%
transmute(CHR, SNP=rsid, BP, A1=effect_allele, A2=neffect_allele, FRQ_A_661866=effect_allele_frequency, FRQ_U_661866=effect_allele_frequency, INFO=1, OR=exp(beta), SE=se, P=p_wald, Direction=effects, HetISqt=0.5, HetDf=1, HetPVa=0.5, Nca=661866, Nco=661866, Neff_half=661866)

daner_stouffer <- 
mc %>%
separate(`chr:pos`, c('CHR', 'BP')) %>%
transmute(CHR, SNP=rsid, BP, A1=effect_allele, A2=neffect_allele, FRQ_A_661866=effect_allele_frequency, FRQ_U_661866=effect_allele_frequency, INFO=1, OR=exp(beta), SE=se, P=p_stouffer, Direction=effects, HetISqt=0.5, HetDf=1, HetPVa=0.5, Nca=661866, Nco=661866, Neff_half=661866)

write_tsv(daner_wald, 'results/meta/metacarpa/pgc_mdd_full_eur_hg19_v3.49.24.mc.wald.gz')
write_tsv(daner_stouffer, 'results/meta/metacarpa/pgc_mdd_full_eur_hg19_v3.49.24.mc.stouffer.gz')


# lambda_gc <- median(qchisq(P, 1, lower.tail=F)) / qchisq(0.5,1)
