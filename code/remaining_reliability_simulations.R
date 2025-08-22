###########################################################
# Finding simulations IDs that did not finish             #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on August 7, 2025                          #
###########################################################

library(tidyverse)

# importing original splits
splits <- read.csv(file = "data/uniform_jaccard_data_splits_1000pairs.csv")

# sanity check
range(splits$split_id)

# listing results files
files <- list.files(
    path = "reliability_jaccard_cluster_results",
    pattern = ".rds", 
    full.names = TRUE
    )

# sorting files numerically by the number in the filename
files <- files[order(readr::parse_number(files) )]

# importing these files
reliability_results <- map2_dfr(
    .x = files,
    .y = seq_along(files),
    .f = ~ readRDS(.x)
    )

# finding splits IDs that are in splits$split_id but not in reliability_results$split_id
splits_in_results <- reliability_results %>%
    mutate(split_id = str_sub(string = split_id, start = 7, end = 11) ) %>%
    pull(split_id) %>%
    unique()
    
all_splits <- sprintf("%05d", splits$split_id) %>% unique()

# splits IDs that are both in splits$split_id and reliability_results$split_id
intersect(splits_in_results, all_splits) %>% length()
intersect(splits_in_results, all_splits)

# splits IDs that are in splits$split_id but not in reliability_results$split_id
setdiff(all_splits, splits_in_results) %>% length()
setdiff(all_splits, splits_in_results)

# exporting a new csv file with remaining simulation IDs
remaining_splits <- splits %>%
    filter(sprintf("%05d", splits$split_id) %in% setdiff(all_splits, splits_in_results) )

# sanity check
remaining_splits %>% pull(split_id) %>% n_distinct()

# exporting it
write.csv(x = remaining_splits, file = "data/uniform_jaccard_data_splits_1000pairs_remaining2.csv", row.names = FALSE)
