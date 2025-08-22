###########################################################
# Defining subsets of participants with uniform Jaccard   #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on July 30, 2025                           #
###########################################################

# importing R packages
library(tidyverse)
library(glue)

# importing the reshaped MEG decoding results
raw_df <- read.csv(file = "data/decoding_results_reshaped.csv") |>
    # removing the last participant (to get a round number of participants)
    mutate(participant_id = cur_group_id(), .by = participant) |>
    dplyr::filter(participant_id < 33) |>
    select(-participant_id)

# getting unique participants
unique_participants <- unique(raw_df$participant)

# getting the number of participants
n_participants <- as.numeric(length(unique_participants) )

# defining the total number of splits
# n_splits <- 50

# creating a dataframe where each participant appears in every split
# split_df <- crossing(participant = unique_participants, split_number = 1:n_splits)

# assigning alternating groups (A or B) within each split for each participant
# split_df <- split_df |>
#     group_by(split_number) |>
#     # shuffling groups for each split
#     mutate(group = sample(rep(c("A", "B"), length.out = n() ) ) ) |>
#     ungroup()

# computing all possible Jaccard values for overlaps 0 to 16
overlaps <- 0:16
jaccard_values <- overlaps / (n_participants - overlaps)
jaccard_values <- round(jaccard_values, 5) # for numerical clarity

jaccard_table <- data.frame(
    overlap = overlaps,
    jaccard = jaccard_values
    )

# checking the results
print(jaccard_table)

# defining some variables
n_total <- n_participants
subset_size <- n_participants / 2
participants <- as.character(1:n_total)

# storing pairs
jaccard_exact_pairs <- vector("list", length(overlaps) )
names(jaccard_exact_pairs) <- paste0("J_", jaccard_values)

# desired number of pairs
n_pairs <- 1000

# getting 1000 pairs of subsets for each level of Jaccard similarity
for (i in seq_along(overlaps) ) {
    
    ov <- overlaps[i]
    pairs <- list()
    count <- 0
    
    while (count < n_pairs) {
        
        shared <- if (ov > 0) sample(participants, ov) else character(0)
        remaining <- setdiff(participants, shared)
        n_to_add <- subset_size - ov
        
        s1 <- sort(c(shared, sample(remaining, n_to_add) ) )
        s2 <- sort(c(shared, sample(setdiff(remaining, s1), n_to_add) ) )
        
        j <- length(intersect(s1, s2) ) / length(union(s1, s2) )
        
        if (round(j, 5) == jaccard_values[i]) {
            
            pairs[[count + 1]] <- list(s1 = s1, s2 = s2, jaccard = j)
            count <- count + 1
            
        }
        
    }
    
    jaccard_exact_pairs[[i]] <- pairs
    
}

# checking the distribution
summary_exact <- data.frame(
    jaccard = names(jaccard_exact_pairs),
    overlap = overlaps,
    n_pairs = sapply(jaccard_exact_pairs, length)
    )

print(summary_exact)

# converting to "split_df" format
# flattening all pairs into one list with associated Jaccard info and split_number
split_list <- list()
split_counter <- 1

for (jaccard_name in names(jaccard_exact_pairs) ) {
    
    for (pair in jaccard_exact_pairs[[jaccard_name]]) {
        
        split_list[[split_counter]] <- tibble(
            participant = c(pair$s1, pair$s2),
            split_id = split_counter,
            group = c(rep("A", length(pair$s1) ), rep("B", length(pair$s2) ) )
            )
        
        split_counter <- split_counter + 1
        
    }
    
}

# combining into a single tibble
split_df <- bind_rows(split_list)

# formatting participant names nicely
split_df <- split_df %>%
    mutate(
        participant = str_pad(participant, 2, pad = "0"),
        participant = paste0("participant", participant)
        ) %>%
    arrange(participant)

# checking the results
# head(split_df, 20)

# final check - computing Jaccard similarity for each split and plotting the distribution
# defining a function to compute Jaccard similarity
# jaccard_sim <- function (a, b) {
#     
#     length(intersect(a, b) ) / length(union(a, b) )
#     
# }

# grouping by split_id and extract A/B groups to compute Jaccard similarity
jaccard_by_split <- split_df %>%
    group_by(split_id) %>%
    summarise(
        group_A = list(participant[group == "A"]),
        group_B = list(participant[group == "B"]),
        # .groups = "drop"
        ) %>%
    # rowwise() %>%
    ungroup() %>%
    # mutate(jaccard = jaccard_sim(group_A[[1]], group_B[[1]]) ) %>%
    # mutate(jaccard = jaccard_sim(group_A[1], group_B[1]) ) %>%
    mutate(
        intersection = map2_int(group_A, group_B, ~length(intersect(.x, .y) ) ),
        union = map2_int(group_A, group_B, ~length(union(.x, .y) ) ),
        jaccard = intersection / union
        )

# sanity checks
table(jaccard_by_split$jaccard)

# plotting the histogram of Jaccard similarities
table(jaccard_by_split$jaccard) %>%
    data.frame() %>%
    ggplot(aes(x = round(as.numeric(Var1), 2), y = Freq) ) +
    geom_bar(stat = "identity") +
    theme_bw() +
    labs(x = "Jaccard similarity", y = "Count")

# creating a function to generate a canonical "signature" for a split: unordered pair {A, B}
create_signature <- function (a, b) {
    
    # converting sorted participant lists to string
    a_str <- paste(sort(a), collapse = "-")
    b_str <- paste(sort(b), collapse = "-")
    
    # sorting the two group strings to make the pair order-invariant
    sorted_pair <- sort(c(a_str, b_str) )
    
    # combining into a single signature
    paste(sorted_pair, collapse = "||")
    
}

# applying signature creation to each split
jaccard_by_split <- jaccard_by_split %>%
    mutate(signature = map2_chr(group_A, group_B, create_signature) )

# checking for potential duplicates
jaccard_by_split %>% nrow()
jaccard_by_split %>% distinct(signature) %>% nrow()

# converting back to the original format of split_df
# split_df2 <- jaccard_by_split %>%
#     select(split_id, group_A, group_B) %>%
#     pivot_longer(group_A:group_B)

# creating two tibbles: one for group A, one for group B
group_A_df <- jaccard_by_split %>%
    select(split_id, group_A) %>%
    unnest(group_A) %>%
    rename(participant = group_A) %>%
    mutate(group = "A")

group_B_df <- jaccard_by_split %>%
    select(split_id, group_B) %>%
    unnest(group_B) %>%
    rename(participant = group_B) %>%
    mutate(group = "B")

# combining and arranging
jaccard_summary <- bind_rows(group_A_df, group_B_df) %>%
    arrange(split_id, participant) %>%
    # re-computing similarity metrics
    # group_by(split_number) %>%
    # mutate(
    #     intersection = map2_int(group_A, group_B, ~length(intersect(.x, .y) ) ),
    #     union = map2_int(group_A, group_B, ~length(union(.x, .y) ) ),
    #     jaccard = intersection / union
    #     ) %>%
    # ungroup()
    group_by(split_id) %>%
    summarise(
        group_A = list(participant[group == "A"]),
        group_B = list(participant[group == "B"]),
        .groups = "drop"
        ) %>%
    mutate(
        intersection = map2_int(group_A, group_B, ~length(intersect(.x, .y) ) ),
        union        = map2_int(group_A, group_B, ~length(union(.x, .y) ) ),
        jaccard      = intersection / union
        ) %>%
    select(split_id, intersection, union, jaccard)

# joining back to the original long-format split_df
split_df_augmented <- split_df %>%
    left_join(jaccard_summary, by = "split_id") %>%
    select(split_id, group, participant, intersection, union, jaccard) %>%
    arrange(split_id, group)

# checking results (544k rows: 17000 pairs of subsets x 32 participants)
head(split_df_augmented, 10)
table(split_df_augmented$jaccard)

# saving this list
write.csv(x = split_df_augmented, file = glue("uniform_jaccard_data_splits_{n_pairs}pairs.csv"), row.names = FALSE)

# splitting existing list
some_df <- read.csv(file = "uniform_jaccard_data_splits_1000pairs.csv")
range(some_df$split_id)
