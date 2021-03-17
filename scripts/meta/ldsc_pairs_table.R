library(dplyr)
library(readr)
library(stringr)

# read in each log text file
ldsc_logs <- lapply(unlist(snakemake@input), readLines)

# find the line of the file with the genetic correlation
# results, then read in the following 2 lines
# as a table
parse_ldsc_log <- function(ldsc_txt) {
	# find start of table and pick lines
	summary_results_title_loc <- which(ldsc_txt == "Summary of Genetic Correlation Results")
	summary_results_header_loc <- summary_results_title_loc + 1
	summary_results_body_loc <- summary_results_title_loc + 2
	
	# extract these lines from the log contents
	summary_results_contents <- ldsc_txt[summary_results_header_loc:summary_results_body_loc]
	
	# convert into a table
	return(read_table2(summary_results_contents))
	
}

# parse LDSC results and merge into a single table
ldsc_results <- 
bind_rows(lapply(ldsc_logs, parse_ldsc_log))

# parse cohorts and subcohorts from sumstats filenames

# regular expression to parse. Format:
# daner_mdd_{COHORT}.{ANCESTRY}.hg19.{SUBCOHORT}.aligned.sumstats.gz
cohort_regx <- "daner_mdd_(\\w+)\\.([a-z]+)\\.hg19\\.(\\w+)\\.aligned\\.sumstats\\.gz"

p1_cohorts <- as_tibble(str_match(basename(ldsc_results$p1), cohort_regx)) %>%
rename(sumstats1=V1, cohort1=V2, subcohort1=V4, ancestry1=V3)

p2_cohorts <- as_tibble(str_match(ldsc_results$p2, cohort_regx))  %>%
rename(sumstats2=V1, cohort2=V2, subcohort2=V4, ancestry2=V3)

cohorts <- bind_cols(p1_cohorts, p2_cohorts)

ldsc_cohort_results <- ldsc_results %>%
mutate(sumstats1=basename(p1), sumstats2=basename(p2)) %>%
inner_join(cohorts, by=c('sumstats1', 'sumstats2')) %>%
select(-p1, -p2, -sumstats1, -sumstats2) %>%
select(cohort1, subcohort1, cohort2, subcohort2, ancestry1, ancestry2, everything()) %>%
arrange(cohort1, subcohort1, cohort2, subcohort2)

write_tsv(ldsc_cohort_results, snakemake@output[[1]])