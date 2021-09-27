# Calculate LDSC covariance structure for GenomicSEM

library(readxl)
library(dplyr)
library(stringr)
library(readr)
library(GenomicSEM)


# Trait sumstats
sumstats <- snakemake@input$sumstats

# traitnames
trait_names=snakemake@params$cohorts

# LD score directory
ld=snakemake@input$w_ld_chr

# Prevalences
# Calculate sample prevalences from Ricopili output
basic_xls <- snakemake@input$samples
names(basic_xls) <- trait_names

basic <- lapply(basic_xls, read_excel)

basic_sum <- bind_rows(lapply(basic, function(b) b %>% filter(Dataset == 'SUM')))

sample_prev <- signif(basic_sum$N_cases / (basic_sum$N_cases + basic_sum$N_controls), 4)

# pop prevalence = 15%
pop_prev = rep(0.15, times=length(sample_prev))

covstruct_capture <-
capture.output(cohorts_covstruct <-
    ldsc(traits=sumstats, trait.names=trait_names,
         sample.prev=sample_prev, population.prev=pop_prev,
         ld=ld, wld=ld))
        
         
# parse ldsc() output
# "Mean Chi"
# "Lambda GC"
# "\"Intercept"
# "Total Observed Scale h2"
# "Total Liability Scale h2"

digit_regx <- "-*[[:digit:]]+\\.[[:digit:]]+"

mean_chisq <- as.numeric(str_extract(str_subset(covstruct_capture, "Mean Chi"), digit_regx))

lambda_gc <- as.numeric(str_extract(str_subset(covstruct_capture, "Lambda GC"), digit_regx))

intercept <- as.numeric(str_extract(str_subset(covstruct_capture, "\"Intercept"), digit_regx))

h2_obs_match <- str_match(sapply(str_split(str_subset(covstruct_capture, "Total Observed Scale h2"), ":"), last), "(.+) \\((.+)\\)")
h2_obs <- as.numeric(h2_obs_match[,2])
h2_se_obs <- as.numeric(h2_obs_match[,3])

h2_liab_match <- str_match(sapply(str_split(str_subset(covstruct_capture, "Total Liability Scale h2"), ":"), last), "(.+) \\((.+)\\)")
h2_liab <- as.numeric(h2_liab_match[,2])
h2_se_liab <- as.numeric(h2_liab_match[,3])


# Output
# deparse 
cohorts_covstruct_r <- snakemake@output$covstruct
dput(cohorts_covstruct, cohorts_covstruct_r, control=c('all', 'digits17'))

ldsc_info <-
tibble(pheno=trait_names, ancestries=snakemake@wildcards$ancestries,
       N_cases=basic_sum$N_cases, N_controls=basic_sum$N_controls,
       sample_prev=sample_prev, pop_prev=pop_prev,
       LambdaGC=basic_sum$`LAMBDA-GC`,
       MeanChiSq=mean_chisq, LambdaGCldsc=lambda_gc, ldsc_intercept=intercept,
       h2_obs, h2_se_obs, h2_liab, h2_se_liab)
       
ldsc_info_txt <- snakemake@output$ldsc_table
write_tsv(ldsc_info, ldsc_info_txt)

