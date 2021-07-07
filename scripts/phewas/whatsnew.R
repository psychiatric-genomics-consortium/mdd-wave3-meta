# compare two ukb releases processed by conv.groovy

library(methods)
library(stringr)
library(readr)
library(dplyr)

# find location of this script
all_args <- commandArgs(trailingOnly = FALSE)
file_flag <- "--file="
script_path <- sub(file_flag, "", all_args[grep(file_flag, all_args)])
script_dir <- dirname(script_path)


# pass in directories where earlier and later releases have been
# extracted

args <- commandArgs(TRUE)
earlier_ukb <- args[1]
later_ukb <- args[2]

# field meta data
showcase <- read_csv(file.path(script_dir, 'Data_Dictionary_Showcase.csv'))
categories <- read_tsv(file.path(script_dir, 'category.txt'))

out_file <- file(file.path(later_ukb, 'WHATSNEW'), 'w')

cat('Earlier release: ', basename(normalizePath(earlier_ukb)), '\n', file=out_file)
cat('Later release: ', basename(normalizePath(later_ukb)), '\n\n\n', file=out_file)

# open list of fields 
print('Opening field lists')
earlier_fields_file <- list.files(earlier_ukb, 'ukb[[:digit:]]*.field', full.name=TRUE)

earlier_fields <- read.table(earlier_fields_file, col.names='field')

later_fields_file <- list.files(later_ukb, 'ukb[[:digit:]]*.fields', full.name=TRUE)
later_fields <- read.table(later_fields_file, col.names='field')

new_fields <- setdiff(later_fields$field, earlier_fields$field)

print('Checking for new fields')
if(length(new_fields) == 0) {
        cat('New fields: None\n', file=out_file)
} else {

        cat('New fields:\n', file=out_file)

        new_showcase <- showcase %>% filter(FieldID %in% new_fields) %>%
        select(FieldID, Category, Field) %>%
        left_join(categories, by=c('Category'='category_id'))

        for(category in unique(new_showcase$Category)) {
                category_fields <- new_showcase %>%filter(Category %in% category)
                cat('*', as.character(category_fields$title[1]), '\n', file=out_file)
                for(i in seq.int(nrow(category_fields))) {
                        cat('* *', as.character(category_fields$FieldID[i]), as.character(category_fields$Field[i]), '\n', file=out_file)
                }
        }

}

# compare data in each field

print('Comparing field data') 

earlier_rds_files <- list.files(earlier_ukb, '.rds', full.name=TRUE)

earlier_rds <- data.frame(earlier_file=earlier_rds_files, rds=sapply(earlier_rds_files, basename))


later_rds_files <- list.files(later_ukb, '.rds', full.name=TRUE)

later_rds <- data.frame(later_file=later_rds_files, rds=sapply(later_rds_files, basename))

rds_files <- later_rds %>%
left_join(earlier_rds, by='rds')

for(i in seq.int(nrow(rds_files))) {

        earlier_file <- as.character(rds_files$earlier_file[i])
        later_file <- as.character(rds_files$later_file[i])
        prefix <- paste(head(str_split(as.character(rds_files$rds[i]), '\\.')[[1]], -1), collapse='.') 

        print(prefix)

       later_data <- readRDS(later_file)

       if(!is.na(earlier_file)) {
          earlier_data <- readRDS(earlier_file)
       } else {
          earlier_data <- later_data
          earlier_data[TRUE] <- NA
       }

       later_data_counts <- plyr::adply(later_data, 2, function(x) sum(!is.na(x)))
       earlier_data_counts <- plyr::adply(earlier_data, 2, function(x) sum(!is.na(x)))
       names(earlier_data_counts) <- c('field', 'earlier_count')
       names(later_data_counts) <- c('field', 'later_count')

       rm(later_data, earlier_data)
       gc()

       data_counts <- later_data_counts %>%
       left_join(earlier_data_counts, by='field') %>%
       filter(field != 'f.eid') %>%
       mutate(earlier_count=if_else(is.na(earlier_count),
                                    true=0L, false=earlier_count),
              FieldID=as.integer(sapply(str_split(field, '\\.'), function(x) x[2]))) %>%
       left_join(showcase %>% select(FieldID, Field), by='FieldID')

       updated_fields <- data_counts %>%
        mutate(new_count=later_count - earlier_count) %>%
        filter(new_count!=0) %>%
        mutate(new_count_str=str_pad(as.character(new_count),
                                     width=max(str_length(as.character(new_count))),
                                     side='left'),
              total_count_str=str_pad(as.character(later_count),
                                     width=max(str_length(as.character(later_count))),
                                     side='left'),
              field_id_str=str_pad(as.character(field),
                                     width=max(str_length(as.character(field))),
                                     side='left')
              )

       if(nrow(updated_fields) > 0) {

               cat(prefix, '\n', file=out_file)
               cat('new / total\n', file=out_file)

               for(i in seq.int(nrow(updated_fields))) {
                       cat('  * ', updated_fields$new_count_str[i], ' / ',
                                   updated_fields$total_count_str[i], 
                                   ' (', updated_fields$field_id_str[i], ') ',
                                   as.character(updated_fields$Field[i]), '\n',
                                   sep='', file=out_file)
               }
            

               cat('\n\n', file=out_file)

       }

}

# bulk fields

earlier_bulk_file <- list.files(earlier_ukb, '.bulk', full.name=TRUE)
later_bulk_file <- list.files(later_ukb, '.bulk', full.name=TRUE)

earlier_bulk <- read.table(earlier_bulk_file, col.names=c('eid', 'field'))

later_bulk <- read.table(later_bulk_file, col.names=c('eid', 'field'))

new_bulk_eids_fields <- 
as_data_frame(later_bulk) %>%
anti_join(earlier_bulk, by=c('eid', 'field'))

total_bulk <- as_data_frame(later_bulk) %>%
group_by(field) %>%
summarize(total=n())

new_bulk_fields <- new_bulk_eids_fields %>%
group_by(field) %>%
summarize(n=n()) %>%
mutate(FieldID=sapply(strsplit(as.character(field), split='_'), function(x) as.integer(x[1]))) %>%
left_join(showcase, by='FieldID') %>%
left_join(total_bulk, by='field') %>%
select(FieldID, Field, n, total) %>%
mutate(new_count_str=str_pad(as.character(n),
                                     width=max(str_length(as.character(n))),
                                     side='left'),
              total_count_str=str_pad(as.character(total),
                                     width=max(str_length(as.character(total))),
                                     side='left'),
              field_id_str=str_pad(paste0('f.', as.character(FieldID)),
                                     width=max(str_length(as.character(FieldID))),
                                     side='left')
              )

cat('BULK fields\n', file=out_file)
cat('new/total\n', file=out_file)

for(i in 1:nrow(new_bulk_fields)) {
        cat('  * ', new_bulk_fields$new_count_str[i], ' / ',
                                   new_bulk_fields$total_count_str[i], 
                                   ' (', new_bulk_fields$field_id_str[i], ') ',
                                   as.character(new_bulk_fields$Field[i]), '\n',
                                   sep='', file=out_file)
}

close(out_file)
