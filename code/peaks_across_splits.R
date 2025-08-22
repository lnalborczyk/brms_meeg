###########################################################
# Defining subsets of participants with uniform Jaccard   #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on July 15, 2025                            #
###########################################################

# importing R packages
library(tidyverse)
library(signal)

# for setting up the cluster
library(doParallel)
library(foreach)

# total number of total CPU cores available on my laptop
total_cores <- 12

# registering parallel backend
cl <- makeCluster(total_cores)
registerDoParallel(cl)

# importing home-made helper functions
source("code/utils.R")

# importing the data splits list
# split_df <- read.csv(file = "uniform_jaccard_data_splits.csv")
split_df <- read.csv(file = "uniform_jaccard_data_splits_100pairs.csv")

# importing the reshaped MEG decoding results
raw_df <- read.csv(file = "data/decoding_results_reshaped.csv") |>
    # removing the last participant (to get a round number of participants)
    mutate(participant_id = cur_group_id(), .by = participant) |>
    dplyr::filter(participant_id < 33) |>
    select(-participant_id)

# averaging across participants
summary_df <- raw_df |>
    summarise(
        auc_mean = mean(auc),
        auc_sd = sd(auc),
        auc_se = sd(auc) / sqrt(n() ),
        .by = c(participant, time)
        )

# creating a dataframe with all (split_number, group) combinations to iterate over
split_group_combinations <- crossing(
    split_number = split_df$split_id,
    group = split_df$group
    ) |>
    # creating a custom simulation ID per group
    mutate(split_id = sprintf("s%05d_%s", split_number, group) )

# for each simulation
reliability_results <- foreach(
    # i = seq_len(n_splits), .combine = bind_rows,
    i = 1:nrow(split_group_combinations), .combine = bind_rows,
    .packages = c("dplyr", "tidyr", "signal")
    ) %dopar% {
    
        # extracting split_number and group for this simulation
        current_split <- split_group_combinations$split_number[i]
        current_group <- split_group_combinations$group[i]
        current_split_id <- split_group_combinations$split_id[i]
        
        # printing progress
        message_parallel(sprintf("\nRunning simulation: %s (%d out of %d)...", current_split_id, i, nrow(split_group_combinations) ) )
        
        # selecting participants from the correct group
        participants_split <- split_df |>
            dplyr::filter(split_id == current_split, group == current_group) |>
            pull(participant)
        
        # subsetting data accordingly
        raw_df_split <- raw_df |> dplyr::filter(participant %in% participants_split)
        summary_df_split <- summary_df |> dplyr::filter(participant %in% participants_split)
        
        # averaging and finding the peak at the group level
        average_meg <- summary_df_split |>
            group_by(time) |>
            summarise(auc_mean = mean(auc_mean) ) |>
            ungroup() |>
            pull(auc_mean)
        
        # smoothing the signal
        smoothed <- sgolayfilt(x = average_meg, p = 3, n = 11)
        
        # finding the the peak
        time <- unique(summary_df_split$time)
        peak_index <- which.max(smoothed)
        peak_time <- time[peak_index]
        
        # plotting
        plot(
            x = time, y = average_meg,
            type = "l", col = "gray",
            # main = "Group-level average AUC with smoothed peak",
            xlab = "Time (s)",
            ylab = "AUC"
            )
        
        lines(time, smoothed, col = "steelblue")
        points(peak_time, smoothed[peak_index], col = "orangered", pch = 19)
        
        # binding all results together
        temp_reliability_results <- data.frame(
            split_id = current_split_id,
            peak_estimate = peak_time
            )
        
        # returning the results
        temp_reliability_results
        
}

# saving the results
write.csv(x = reliability_results, file = "peaks_across_splits_100pairs.csv", row.names = FALSE)

# stopping the cluster
stopCluster(cl)
