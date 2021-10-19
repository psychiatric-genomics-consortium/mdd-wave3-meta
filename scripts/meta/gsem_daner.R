library(readr)
library(dplyr)

gsem <- read_tsv(snakemake@input$gsem)
daner <- read_tsv(snakemake@input$daner)

# join with daner data and turn alleles to match daner
gsem_daner <- gsem %>%
inner_join(daner, by=c('CHR', 'SNP', 'BP'), suffix=c('.gsem', '.daner')) %>%
filter((A1.gsem == A1.daner & A2.gsem == A2.daner) |
        A1.gsem == A2.daner & A2.gsem == A1.daner) %>%
mutate(BETA=if_else(A1.gsem == A1.daner, true=est, false=-est)) %>%
select(CHR, SNP, BP, A1=A1.daner, A2=A2.daner, starts_with('FRQ_A'), starts_with('FRQ_U'),
       INFO, BETA, SE=se_c, P=Pval_Estimate,
       ngt, Direction,
       HetISqt=Q, HetDf=Q_df, HetPVa=Q_pval,
       Nca, Nco, Neff_half)

write_tsv(gsem_daner, snakemake@output[[1]])