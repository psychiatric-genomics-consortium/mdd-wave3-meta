library(dplyr)
library(readr)
library(tidyr)

carpa_txt <- snakemake@input$metacarpa
carpa <- read_tsv(carpa_txt, col_types=list(`chr:pos`=col_character()))

daner_rp_gz <- snakemake@input$ricopili
daner_rp <- read_tsv(daner_rp_gz, col_types=list("SNP"=col_character()))

# find largest sample size
Nco <- Nca <- round(max(carpa$n)/2)

frq_a_col <- paste('FRQ', 'A', Nca, sep='_')
frq_u_col <- paste('FRQ', 'U', Nca, sep='_')

carpa_rp <- 
carpa %>%
inner_join(daner_rp, by=c('rsid'='SNP'), suffix=c('.mc', '.rp'))

alleles_non_matching <- carpa_rp %>% filter(effect_allele != A1 | neffect_allele != A2)

if(nrow(alleles_non_matching > 0)) {
    cat("Non-matching effect/noneffect alleles when merging Metacarpa and Ricopili results\n")
    stop()
}

daner_mc <-
carpa_rp %>%
mutate(or_wald=exp(beta)) %>%
select(CHR, SNP=rsid, BP, A1, A2,
    starts_with('FRQ_A'), starts_with('FRQ_U'), INFO,
    OR=or_wald, SE=se, P=p_wald, ngt,
    Direction=effects, HetISqt, HetDf, HetPVa, 
    Nca, Nco, Neff_half)

daner_gz <- snakemake@output[[1]]
write_tsv(daner_mc, daner_gz)

# lambda_gc <- median(qchisq(P, 1, lower.tail=F)) / qchisq(0.5,1)
