# merge LDSC log tables with Open GWAS info

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# Read in LDSC rg tables

ldsc_cols <- cols(
  p1 = col_character(),
  p2 = col_character(),
  rg = col_double(),
  se = col_double(),
  z = col_double(),
  p = col_double(),
  h2_obs = col_double(),
  h2_obs_se = col_double(),
  h2_int = col_double(),
  h2_int_se = col_double(),
  gcov_int = col_double(),
  gcov_int_se = col_double()
)

ldsc_rg <- bind_rows(lapply(snakemake@input$ldsc_full, read_table2, col_types=ldsc_cols))
ldsc_noukbb_rg <- bind_rows(lapply(snakemake@input$ldsc_noukbb, read_table2, col_types=ldsc_cols))
ldsc_previous_rg <-  bind_rows(lapply(snakemake@input$ldsc_previous, read_table2, col_types=ldsc_cols))


# Read in GWAS info

gwasinfo <- read_tsv(snakemake@input$gwasinfo)

# Remove results with `NA`s and merge with GWAS info
ldsc_rg_info <- bind_rows(ldsc_rg, ldsc_noukbb_rg, ldsc_previous_rg) %>%
mutate(id=sapply(str_split(basename(p2), "\\."), first)) %>%
left_join(gwasinfo, by='id') %>%
filter(!is.na(p))


# FDR corrected associations for full results, removing phenotypes from UKB and with large genetic covariance intercepts
ldsc_rg_mr_candidates <- 
ldsc_rg_info %>%
filter(str_detect(p1, pattern='pgc_mdd_full')) %>%
filter(abs(gcov_int) <= 0.05, !str_detect(id, 'ukb')) %>%
mutate(qvalue=fdrtool::fdrtool(p, statistic='p', plot=FALSE)$qval) %>%
filter(qvalue <= 0.05) %>%
arrange(desc(abs(rg))) %>%
select(id, trait, rg, p, qvalue, gcov_int, subcategory)

write_tsv(ldsc_rg_info, 'docs/tables/ldsc_open_rg.txt')
write_tsv(ldsc_rg_mr_candidates, 'docs/tables/ldsc_open_mr_candidates.txt')




