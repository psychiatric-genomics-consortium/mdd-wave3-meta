library(readr)
library(dplyr)

auto_clump <- read_table(snakemake@input$auto,
    col_types = cols("ngt" = col_character()))
allo_clump <- read_table(snakemake@input$allo,
    col_types = cols("ngt" = col_character()))

# get names of FRQ_A and FRQ_U columns
frq_a_auto_col <- names(select(auto_clump, starts_with('FRQ_A')))
frq_u_auto_col <- names(select(auto_clump, starts_with('FRQ_U')))
frq_a_allo_col <- names(select(allo_clump, starts_with('FRQ_A')))
frq_u_allo_col <- names(select(allo_clump, starts_with('FRQ_U')))

clump <- bind_rows(
auto_clump,
allo_clump %>% rename(!!frq_a_auto_col:=frq_a_allo_col,
                      !!frq_u_auto_col:=frq_u_allo_col)
) %>%
mutate(ngt = if_else(ngt == '-', true = 0, false = as.numeric(ngt))) %>%
filter(P <= 5e-8) %>%
arrange(P)

write_tsv(clump, snakemake@output[[1]])