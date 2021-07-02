library(dplyr)
library(readr)
library(tidyr)
library(stringr)

metal_tbl <- snakemake@input$metal
metal <- read_tsv(metal_tbl, col_types=list(MarkerName=col_character()))

daner_rp_gz <- snakemake@input$ricopili
daner_rp <- read_tsv(daner_rp_gz, col_types=list("SNP"=col_character()))


# OR = exp(Z / sqrt(Neff/4*FRQ*(1-FRQ)))
# SE = 1/sqrt(Neff/4*FRQ*(1-MAF))	

metal_rp <- 
metal %>%
inner_join(daner_rp, by=c('MarkerName'='SNP'), suffix=c('.ml', '.rp'))

alleles_non_matching <- 
metal_rp %>%
    filter(!(str_to_upper(Allele1) == A1 & str_to_upper(Allele2) == A2) &
           !(str_to_upper(Allele1) == A2 & str_to_upper(Allele2) == A1))  

if(nrow(alleles_non_matching > 0)) {
    cat("Non-matching alleles when merging Metal and Ricopili results\n")
    stop()
}

# get names of FRQ_A and FRQ_U columns
frq_a_col <- names(select(metal_rp, starts_with('FRQ_A')))
frq_u_col <- names(select(metal_rp, starts_with('FRQ_U')))

daner_ml <-
metal_rp %>%
filter(N > 1) %>%
mutate(ALLELE1=str_to_upper(Allele1),
       ALLELE2=str_to_upper(Allele2),
       frq=.data[[frq_a_col]]) %>%
mutate(Z=case_when(ALLELE1 == A1 & ALLELE2 == A2 ~ Zscore,
                   ALLELE1 == A2 & ALLELE2 == A1 ~ -Zscore,
                   TRUE ~ NA_real_)) %>%
mutate(or_metal=exp(Z/sqrt(N/4*frq*(1-frq))),
       se_metal=1/sqrt(N/4*frq*(1-frq)),
       N_half=N/2) %>%
select(CHR, SNP=MarkerName, BP, A1, A2,
    starts_with('FRQ_A'), starts_with('FRQ_U'), INFO,
    OR=or_metal, SE=se_metal, P=`P-value`, Direction=Direction.rp,
    HetISqt=HetISq, HetDf=HetDf.ml, HetPVa=HetPVal, 
    Nca, Nco, Neff_half=N_half) %>%
arrange(CHR, BP)

daner_gz <- snakemake@output[[1]]
write_tsv(daner_ml, daner_gz)

# lambda_gc <- median(qchisq(P, 1, lower.tail=F)) / qchisq(0.5,1)
